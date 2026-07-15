#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day020_confined_water_axial_radial_density/"
    "regional_classification"
)

INPUT_CSV = (
    INPUT_ROOT
    / "confined_water_region_sensitivity.csv"
)

OUTPUT_CSV = (
    INPUT_ROOT
    / "confined_water_region_sensitivity_audit.csv"
)

OUTPUT_MD = (
    INPUT_ROOT
    / "CONFINED_WATER_REGION_SENSITIVITY_AUDIT_DAY020.md"
)

REGION_ORDER = (
    "interior_core",
    "interfacial_shell",
    "mouth_zones",
    "exterior",
)

DISPLAY_NAMES = {
    "interior_core": "Interior core",
    "interfacial_shell": "Interfacial shell",
    "mouth_zones": "Mouth zones",
    "exterior": "Exterior",
}

CENTRAL_INTERFACE_WIDTH_NM = 0.25
CENTRAL_MOUTH_WIDTH_NM = 0.50


def read_rows() -> list[dict[str, str]]:
    if not INPUT_CSV.exists():
        raise RuntimeError(
            f"Missing sensitivity input: {INPUT_CSV}"
        )

    with INPUT_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(csv.DictReader(handle))


def write_csv(
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            "No sensitivity-audit rows generated"
        )

    with OUTPUT_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rows = read_rows()

    grouped: dict[
        str,
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in rows:
        grouped[row["region"]].append(row)

    missing = [
        region
        for region in REGION_ORDER
        if region not in grouped
    ]

    if missing:
        raise RuntimeError(
            "Missing regions: "
            + ", ".join(missing)
        )

    audit_rows: list[
        dict[str, object]
    ] = []

    for region in REGION_ORDER:
        region_rows = grouped[region]

        counts = [
            float(
                row[
                    "average_water_count"
                ]
            )
            for row in region_rows
        ]

        fractions = [
            float(
                row[
                    "water_count_fraction"
                ]
            )
            for row in region_rows
        ]

        densities = [
            float(
                row[
                    "volume_weighted_density_nm^-3"
                ]
            )
            for row in region_rows
        ]

        central_matches = [
            row
            for row in region_rows
            if abs(
                float(
                    row[
                        "interface_half_width_nm"
                    ]
                )
                - CENTRAL_INTERFACE_WIDTH_NM
            ) < 1.0e-12
            and abs(
                float(
                    row[
                        "mouth_half_width_nm"
                    ]
                )
                - CENTRAL_MOUTH_WIDTH_NM
            ) < 1.0e-12
        ]

        if len(central_matches) != 1:
            raise RuntimeError(
                f"Could not uniquely identify "
                f"central condition for {region}"
            )

        central = central_matches[0]

        central_count = float(
            central[
                "average_water_count"
            ]
        )

        count_minimum = min(counts)
        count_maximum = max(counts)

        fraction_minimum = min(
            fractions
        )

        fraction_maximum = max(
            fractions
        )

        count_span = (
            count_maximum
            - count_minimum
        )

        relative_span = (
            count_span
            / central_count
            if central_count > 0.0
            else float("nan")
        )

        audit_rows.append(
            {
                "region": region,
                "display_name": (
                    DISPLAY_NAMES[
                        region
                    ]
                ),
                "n_sensitivity_conditions": (
                    len(region_rows)
                ),
                "central_count": (
                    central_count
                ),
                "minimum_count": (
                    count_minimum
                ),
                "maximum_count": (
                    count_maximum
                ),
                "absolute_count_span": (
                    count_span
                ),
                "relative_count_span": (
                    relative_span
                ),
                "central_fraction": float(
                    central[
                        "water_count_fraction"
                    ]
                ),
                "minimum_fraction": (
                    fraction_minimum
                ),
                "maximum_fraction": (
                    fraction_maximum
                ),
                "fraction_span": (
                    fraction_maximum
                    - fraction_minimum
                ),
                "minimum_density_nm^-3": (
                    min(densities)
                ),
                "maximum_density_nm^-3": (
                    max(densities)
                ),
            }
        )

    write_csv(
        audit_rows
    )

    with OUTPUT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Confined-Water "
            "Region Sensitivity Audit\n\n"
        )

        handle.write(
            "The regional partition was evaluated "
            "for nine combinations of interfacial "
            "half-width and mouth half-width.\n\n"
        )

        handle.write(
            "| Region | Central count | "
            "Count range | Fraction range | "
            "Relative span |\n"
        )

        handle.write(
            "|---|---:|---:|---:|---:|\n"
        )

        for row in audit_rows:
            handle.write(
                f"| {row['display_name']} "
                f"| {row['central_count']:.3f} "
                f"| {row['minimum_count']:.3f}–"
                f"{row['maximum_count']:.3f} "
                f"| {100.0 * row['minimum_fraction']:.3f}–"
                f"{100.0 * row['maximum_fraction']:.3f}% "
                f"| {100.0 * row['relative_count_span']:.3f}% |\n"
            )

        handle.write(
            "\n## Interpretation\n\n"
        )

        handle.write(
            "Large variation in a regional count "
            "indicates that the result depends strongly "
            "on the operational boundary. Such a region "
            "should not be reported as a unique physical "
            "population until the density profiles are "
            "used to identify structural minima.\n"
        )

    print(
        "Day020 confined-water regional "
        "sensitivity audit completed.",
        flush=True,
    )

    print(
        f"Input rows: {len(rows)}",
        flush=True,
    )

    print(
        f"Regions: {len(audit_rows)}",
        flush=True,
    )

    for row in audit_rows:
        print(
            f"{row['display_name']}: "
            f"{row['minimum_count']:.3f} to "
            f"{row['maximum_count']:.3f} waters; "
            f"relative span "
            f"{100.0 * row['relative_count_span']:.2f}%",
            flush=True,
        )

    print(
        "Overall validation: PASS",
        flush=True,
    )

    print(
        "Wrote: "
        f"{INPUT_ROOT.relative_to(PROJECT_ROOT)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
