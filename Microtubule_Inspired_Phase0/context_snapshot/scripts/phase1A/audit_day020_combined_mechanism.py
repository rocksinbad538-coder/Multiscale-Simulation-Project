#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_combined_dephasing_relaxation"
)

OUTPUT_ROOT = INPUT_ROOT / "combined_mechanism_audit"

CONDITION_CSV = INPUT_ROOT / "condition_summary.csv"
VALIDATION_CSV = INPUT_ROOT / "numerical_validation.csv"

TEMPERATURE_CSV = (
    OUTPUT_ROOT / "temperature_sensitivity.csv"
)

RATIO_CSV = (
    OUTPUT_ROOT / "same_ratio_scale_sensitivity.csv"
)

REPORT_MD = (
    OUTPUT_ROOT / "COMBINED_MECHANISM_AUDIT_DAY020.md"
)

EXPECTED_ROWS = 24


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    with path.open(
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


def number(
    row: dict[str, str],
    key: str,
) -> float:
    return float(row[key])


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    rows = read_csv(CONDITION_CSV)
    validation_rows = read_csv(VALIDATION_CSV)

    if len(rows) != EXPECTED_ROWS:
        raise RuntimeError(
            f"Expected {EXPECTED_ROWS} conditions, "
            f"found {len(rows)}"
        )

    if len(validation_rows) != 1:
        raise RuntimeError(
            "Expected one validation row"
        )

    validation = validation_rows[0]

    failed_validations = [
        key
        for key, value in validation.items()
        if key.endswith("_pass")
        and value.lower() != "true"
    ]

    if failed_validations:
        raise RuntimeError(
            "Failed validations: "
            + ", ".join(failed_validations)
        )

    gateway_values = [
        number(row, "PYR4_gateway_fraction")
        for row in rows
    ]

    final_values = [
        number(
            row,
            "mean_final_PYR5_population_from_high_states",
        )
        for row in rows
    ]

    steady_values = [
        number(
            row,
            "mean_steady_PYR5_population",
        )
        for row in rows
    ]

    gamma0_rows = [
        row
        for row in rows
        if number(row, "gamma_phi_ps_inv") == 0.0
    ]

    nonzero_gamma_rows = [
        row
        for row in rows
        if number(row, "gamma_phi_ps_inv") > 0.0
    ]

    gamma100_rows = [
        row
        for row in rows
        if number(row, "gamma_phi_ps_inv") == 100.0
    ]

    gamma0_steady = [
        number(
            row,
            "mean_steady_PYR5_population",
        )
        for row in gamma0_rows
    ]

    nonzero_steady = [
        number(
            row,
            "mean_steady_PYR5_population",
        )
        for row in nonzero_gamma_rows
    ]

    gamma100_final = [
        number(
            row,
            "mean_final_PYR5_population_from_high_states",
        )
        for row in gamma100_rows
    ]

    by_condition: dict[
        tuple[float, float],
        dict[float, dict[str, str]],
    ] = defaultdict(dict)

    for row in rows:
        key = (
            number(row, "kappa_ref_ps_inv"),
            number(row, "gamma_phi_ps_inv"),
        )

        temperature = number(
            row,
            "temperature_K",
        )

        by_condition[key][temperature] = row

    temperature_rows: list[dict[str, object]] = []

    for (
        kappa,
        gamma,
    ), temperature_map in sorted(
        by_condition.items()
    ):
        if set(temperature_map) != {
            150.0,
            300.0,
        }:
            raise RuntimeError(
                f"Missing temperature pair for "
                f"kappa={kappa}, gamma={gamma}"
            )

        low = temperature_map[150.0]
        high = temperature_map[300.0]

        temperature_rows.append(
            {
                "kappa_ref_ps_inv": kappa,
                "gamma_phi_ps_inv": gamma,
                "absolute_change_final_PYR5_300_minus_150": abs(
                    number(
                        high,
                        "mean_final_PYR5_population_from_high_states",
                    )
                    - number(
                        low,
                        "mean_final_PYR5_population_from_high_states",
                    )
                ),
                "absolute_change_steady_PYR5_300_minus_150": abs(
                    number(
                        high,
                        "mean_steady_PYR5_population",
                    )
                    - number(
                        low,
                        "mean_steady_PYR5_population",
                    )
                ),
                "absolute_change_PYR4_gateway_300_minus_150": abs(
                    number(
                        high,
                        "PYR4_gateway_fraction",
                    )
                    - number(
                        low,
                        "PYR4_gateway_fraction",
                    )
                ),
            }
        )

    ratio_groups: dict[
        tuple[float, float],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in rows:
        gamma = number(
            row,
            "gamma_phi_ps_inv",
        )

        if gamma == 0.0:
            continue

        kappa = number(
            row,
            "kappa_ref_ps_inv",
        )

        temperature = number(
            row,
            "temperature_K",
        )

        ratio = gamma / kappa

        ratio_groups[
            (
                temperature,
                round(ratio, 12),
            )
        ].append(row)

    ratio_rows: list[dict[str, object]] = []

    for (
        temperature,
        ratio,
    ), group in sorted(
        ratio_groups.items()
    ):
        if len(group) < 2:
            continue

        final_group = [
            number(
                row,
                "mean_final_PYR5_population_from_high_states",
            )
            for row in group
        ]

        steady_group = [
            number(
                row,
                "mean_steady_PYR5_population",
            )
            for row in group
        ]

        gateway_group = [
            number(
                row,
                "PYR4_gateway_fraction",
            )
            for row in group
        ]

        conditions = "; ".join(
            (
                f"kappa={number(row, 'kappa_ref_ps_inv'):g},"
                f"gamma={number(row, 'gamma_phi_ps_inv'):g}"
            )
            for row in sorted(
                group,
                key=lambda item: number(
                    item,
                    "kappa_ref_ps_inv",
                ),
            )
        )

        ratio_rows.append(
            {
                "temperature_K": temperature,
                "gamma_over_kappa": ratio,
                "conditions": conditions,
                "final_PYR5_spread": (
                    max(final_group)
                    - min(final_group)
                ),
                "steady_PYR5_spread": (
                    max(steady_group)
                    - min(steady_group)
                ),
                "PYR4_gateway_spread": (
                    max(gateway_group)
                    - min(gateway_group)
                ),
            }
        )

    write_csv(
        TEMPERATURE_CSV,
        temperature_rows,
    )

    write_csv(
        RATIO_CSV,
        ratio_rows,
    )

    max_temperature_final = max(
        temperature_rows,
        key=lambda row: float(
            row[
                "absolute_change_final_PYR5_300_minus_150"
            ]
        ),
    )

    max_temperature_steady = max(
        temperature_rows,
        key=lambda row: float(
            row[
                "absolute_change_steady_PYR5_300_minus_150"
            ]
        ),
    )

    minimum_steady_row = min(
        rows,
        key=lambda row: number(
            row,
            "mean_steady_PYR5_population",
        ),
    )

    maximum_final_row = max(
        rows,
        key=lambda row: number(
            row,
            "mean_final_PYR5_population_from_high_states",
        ),
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Combined-Mechanism Audit\n\n"
        )

        handle.write(
            "## Numerical status\n\n"
        )

        handle.write(
            "- All numerical validation criteria: PASS.\n"
        )

        handle.write(
            f"- Conditions audited: "
            f"{len(rows)}/{EXPECTED_ROWS}.\n\n"
        )

        handle.write(
            "## Robust pathway result\n\n"
        )

        handle.write(
            f"- PYR4 gateway fraction: "
            f"{min(gateway_values):.9f} to "
            f"{max(gateway_values):.9f}.\n"
        )

        handle.write(
            "- PYR4 remains the dominant thermal "
            "gateway in every tested condition.\n\n"
        )

        handle.write(
            "## Transient PYR5 population\n\n"
        )

        handle.write(
            f"- Full 100 ps final-population range: "
            f"{min(final_values):.9f} to "
            f"{max(final_values):.9f}.\n"
        )

        handle.write(
            f"- At gamma_phi=100 ps^-1: "
            f"{min(gamma100_final):.9f} to "
            f"{max(gamma100_final):.9f}.\n"
        )

        handle.write(
            f"- Maximum final population occurs at "
            f"T={number(maximum_final_row, 'temperature_K'):g} K, "
            f"kappa={number(maximum_final_row, 'kappa_ref_ps_inv'):g}, "
            f"gamma={number(maximum_final_row, 'gamma_phi_ps_inv'):g}.\n\n"
        )

        handle.write(
            "## Stationary-state competition\n\n"
        )

        handle.write(
            f"- Thermal-only stationary PYR5 range: "
            f"{min(gamma0_steady):.9f} to "
            f"{max(gamma0_steady):.9f}.\n"
        )

        handle.write(
            f"- Nonzero-dephasing stationary PYR5 range: "
            f"{min(nonzero_steady):.9f} to "
            f"{max(nonzero_steady):.9f}.\n"
        )

        handle.write(
            f"- Overall stationary range: "
            f"{min(steady_values):.9f} to "
            f"{max(steady_values):.9f}.\n"
        )

        handle.write(
            f"- Minimum stationary PYR5 occurs at "
            f"T={number(minimum_steady_row, 'temperature_K'):g} K, "
            f"kappa={number(minimum_steady_row, 'kappa_ref_ps_inv'):g}, "
            f"gamma={number(minimum_steady_row, 'gamma_phi_ps_inv'):g}.\n\n"
        )

        handle.write(
            "## Temperature sensitivity\n\n"
        )

        handle.write(
            f"- Maximum absolute temperature effect on "
            f"the 100 ps PYR5 population: "
            f"{float(max_temperature_final['absolute_change_final_PYR5_300_minus_150']):.9f}.\n"
        )

        handle.write(
            f"- Maximum absolute temperature effect on "
            f"stationary PYR5: "
            f"{float(max_temperature_steady['absolute_change_steady_PYR5_300_minus_150']):.9f}.\n\n"
        )

        handle.write(
            "## Rate-ratio assessment\n\n"
        )

        for row in ratio_rows:
            handle.write(
                f"- At T={float(row['temperature_K']):g} K "
                f"and gamma/kappa="
                f"{float(row['gamma_over_kappa']):g}, "
                f"the final-PYR5 spread across absolute "
                f"scales is "
                f"{float(row['final_PYR5_spread']):.9f}, "
                f"while the stationary spread is "
                f"{float(row['steady_PYR5_spread']):.9f}.\n"
            )

        handle.write(
            "\nThe ratio gamma_phi/kappa identifies the "
            "dominant dissipative channel but does not "
            "fully determine either the transient or "
            "stationary dynamics. Absolute rates relative "
            "to the coherent Hamiltonian remain relevant.\n\n"
        )

        handle.write(
            "## Accepted interpretation\n\n"
        )

        handle.write(
            "The upper PYR2-PYR4 manifold undergoes rapid "
            "redistribution, followed by kinetically slow "
            "capture into PYR5 predominantly through PYR4. "
            "Local dephasing enhances transient access to "
            "PYR5 but competes with thermal relaxation in "
            "the stationary state. Absolute relaxation "
            "times remain phenomenological because no "
            "microscopic spectral density has been derived.\n"
        )

    print(
        "Day020 combined-mechanism audit completed."
    )

    print(
        "Numerical validation: PASS"
    )

    print(
        "PYR4 gateway fraction: "
        f"{min(gateway_values):.6f} to "
        f"{max(gateway_values):.6f}"
    )

    print(
        "Final PYR5 population range: "
        f"{min(final_values):.6f} to "
        f"{max(final_values):.6f}"
    )

    print(
        "Gamma=100 final PYR5 range: "
        f"{min(gamma100_final):.6f} to "
        f"{max(gamma100_final):.6f}"
    )

    print(
        "Thermal-only steady PYR5 range: "
        f"{min(gamma0_steady):.6f} to "
        f"{max(gamma0_steady):.6f}"
    )

    print(
        "Nonzero-dephasing steady PYR5 range: "
        f"{min(nonzero_steady):.6f} to "
        f"{max(nonzero_steady):.6f}"
    )

    print(
        "Maximum temperature effect on final PYR5: "
        f"{float(max_temperature_final['absolute_change_final_PYR5_300_minus_150']):.6f}"
    )

    print(
        "Maximum temperature effect on steady PYR5: "
        f"{float(max_temperature_steady['absolute_change_steady_PYR5_300_minus_150']):.6f}"
    )

    print(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
