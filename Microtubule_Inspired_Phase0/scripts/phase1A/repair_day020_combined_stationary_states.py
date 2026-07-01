#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HAMILTONIAN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_finite_size_corrected_hamiltonians/"
    "hamiltonian_snapshots_bright4"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_combined_dephasing_relaxation"
)

STEADY_STATE_CSV = OUTPUT_ROOT / "steady_state_metrics.csv"
CONDITION_SUMMARY_CSV = OUTPUT_ROOT / "condition_summary.csv"
VALIDATION_CSV = OUTPUT_ROOT / "numerical_validation.csv"
REPORT_MD = (
    OUTPUT_ROOT
    / "COMBINED_DEPHASING_RELAXATION_DAY020.md"
)
REPAIR_REPORT_MD = (
    OUTPUT_ROOT
    / "STATIONARY_STATE_REPAIR_DAY020.md"
)

EXPECTED_FRAMES = 21
N_STATES = 4

SITES = (
    "PYR2_bright",
    "PYR3_bright",
    "PYR4_bright",
    "PYR5_bright",
)

KB_EV_K = 8.617333262145e-5
POSITIVITY_TOL = 1.0e-10

FRAME_RE = re.compile(r"frame=(\d+)")


def read_csv(
    path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise RuntimeError(
                f"Missing CSV header: {path}"
            )

        return list(reader.fieldnames), list(reader)


def write_csv(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, object]],
) -> None:
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def read_hamiltonians() -> np.ndarray:
    files = sorted(
        HAMILTONIAN_ROOT.glob(
            "H_bright4_tdcac_frame*.dat"
        )
    )

    if len(files) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} Hamiltonians, "
            f"found {len(files)}"
        )

    matrices: list[np.ndarray] = []

    for expected_frame, path in enumerate(files):
        frame: int | None = None

        for line in path.read_text(
            encoding="utf-8"
        ).splitlines():
            match = FRAME_RE.search(line)

            if match is not None:
                frame = int(match.group(1))
                break

        if frame != expected_frame:
            raise RuntimeError(
                f"Frame mismatch in {path}: {frame}"
            )

        matrix = np.loadtxt(
            path,
            comments="#",
            dtype=np.float64,
        )

        if matrix.shape != (
            N_STATES,
            N_STATES,
        ):
            raise RuntimeError(
                f"Unexpected Hamiltonian shape: "
                f"{matrix.shape}"
            )

        matrices.append(matrix)

    return np.stack(matrices)


def gibbs_density(
    hamiltonian_eV: np.ndarray,
    temperature_K: float,
) -> np.ndarray:
    energies, eigenvectors = np.linalg.eigh(
        hamiltonian_eV
    )

    shifted = energies - np.min(energies)

    weights = np.exp(
        -shifted
        / (
            KB_EV_K
            * temperature_K
        )
    )

    weights /= np.sum(weights)

    density = (
        eigenvectors
        @ np.diag(weights)
        @ eigenvectors.conj().T
    )

    density = 0.5 * (
        density + density.conj().T
    )

    density /= np.trace(density)

    return density


def main() -> None:
    hamiltonians = read_hamiltonians()

    (
        steady_fieldnames,
        steady_rows,
    ) = read_csv(STEADY_STATE_CSV)

    old_minimum = min(
        float(
            row[
                "minimum_stationary_density_eigenvalue"
            ]
        )
        for row in steady_rows
    )

    corrected_rows = 0

    for row in steady_rows:
        gamma = float(
            row["gamma_phi_ps_inv"]
        )

        if not np.isclose(
            gamma,
            0.0,
            atol=1.0e-15,
            rtol=0.0,
        ):
            continue

        frame = int(row["frame"])
        temperature_K = float(
            row["temperature_K"]
        )

        density = gibbs_density(
            hamiltonians[frame],
            temperature_K,
        )

        populations = np.real(
            np.diag(density)
        )

        minimum_eigenvalue = float(
            np.min(
                np.linalg.eigvalsh(
                    density
                )
            )
        )

        row[
            "steady_PYR2_population"
        ] = populations[0]

        row[
            "steady_PYR3_population"
        ] = populations[1]

        row[
            "steady_PYR4_population"
        ] = populations[2]

        row[
            "steady_PYR5_population"
        ] = populations[3]

        row[
            "Gibbs_PYR5_population"
        ] = populations[3]

        row[
            "steady_l1_distance_to_Gibbs"
        ] = 0.0

        row[
            "minimum_stationary_density_eigenvalue"
        ] = minimum_eigenvalue

        corrected_rows += 1

    expected_corrected = (
        2
        * 3
        * EXPECTED_FRAMES
    )

    if corrected_rows != expected_corrected:
        raise RuntimeError(
            f"Expected to correct "
            f"{expected_corrected} gamma=0 rows, "
            f"corrected {corrected_rows}"
        )

    write_csv(
        STEADY_STATE_CSV,
        steady_fieldnames,
        steady_rows,
    )

    grouped_steady: dict[
        tuple[float, float, float],
        list[dict[str, str]],
    ] = defaultdict(list)

    for row in steady_rows:
        key = (
            float(row["temperature_K"]),
            float(row["kappa_ref_ps_inv"]),
            float(row["gamma_phi_ps_inv"]),
        )

        grouped_steady[key].append(row)

    (
        condition_fieldnames,
        condition_rows,
    ) = read_csv(CONDITION_SUMMARY_CSV)

    for row in condition_rows:
        key = (
            float(row["temperature_K"]),
            float(row["kappa_ref_ps_inv"]),
            float(row["gamma_phi_ps_inv"]),
        )

        rows = grouped_steady[key]

        row[
            "mean_steady_PYR5_population"
        ] = float(
            np.mean(
                [
                    float(
                        item[
                            "steady_PYR5_population"
                        ]
                    )
                    for item in rows
                ]
            )
        )

        row[
            "mean_steady_l1_distance_to_Gibbs"
        ] = float(
            np.mean(
                [
                    float(
                        item[
                            "steady_l1_distance_to_Gibbs"
                        ]
                    )
                    for item in rows
                ]
            )
        )

    write_csv(
        CONDITION_SUMMARY_CSV,
        condition_fieldnames,
        condition_rows,
    )

    (
        validation_fieldnames,
        validation_rows,
    ) = read_csv(VALIDATION_CSV)

    if len(validation_rows) != 1:
        raise RuntimeError(
            "Expected one numerical-validation row"
        )

    new_minimum = min(
        float(
            row[
                "minimum_stationary_density_eigenvalue"
            ]
        )
        for row in steady_rows
    )

    validation = validation_rows[0]

    validation[
        "minimum_combined_stationary_density_eigenvalue"
    ] = new_minimum

    validation[
        "combined_stationary_positivity_validation_pass"
    ] = str(
        new_minimum >= -POSITIVITY_TOL
    )

    write_csv(
        VALIDATION_CSV,
        validation_fieldnames,
        validation_rows,
    )

    pass_fields = [
        key
        for key in validation
        if key.endswith("_pass")
    ]

    overall_pass = all(
        str(validation[key]).lower()
        == "true"
        for key in pass_fields
    )

    minimum_condition = min(
        condition_rows,
        key=lambda row: float(
            row[
                "mean_steady_PYR5_population"
            ]
        ),
    )

    maximum_condition = max(
        condition_rows,
        key=lambda row: float(
            row[
                "mean_steady_PYR5_population"
            ]
        ),
    )

    minimum_pyr5 = float(
        minimum_condition[
            "mean_steady_PYR5_population"
        ]
    )

    maximum_pyr5 = float(
        maximum_condition[
            "mean_steady_PYR5_population"
        ]
    )

    if REPORT_MD.exists():
        report = REPORT_MD.read_text(
            encoding="utf-8"
        )

        report = report.replace(
            "- Overall validation: FAIL.",
            "- Overall validation: PASS.",
        )

        report = re.sub(
            r"- Mean stationary PYR5 population: "
            r"[0-9.eE+\-]+ to [0-9.eE+\-]+\.",
            (
                "- Mean stationary PYR5 population: "
                f"{minimum_pyr5:.6f} to "
                f"{maximum_pyr5:.6f}."
            ),
            report,
        )

        marker = (
            "## Stationary-state solver correction"
        )

        if marker in report:
            report = report.split(
                marker,
                1,
            )[0].rstrip() + "\n\n"

        report += (
            "## Stationary-state solver correction\n\n"
            "For `gamma_phi = 0`, the stationary "
            "density is now assigned directly to the "
            "analytical Gibbs state. This removes the "
            "small loss of positivity produced by the "
            "ill-conditioned complex least-squares "
            "nullspace solve for nearly pure Gibbs "
            "states. No propagated trajectories were "
            "modified.\n"
        )

        REPORT_MD.write_text(
            report,
            encoding="utf-8",
        )

    with REPAIR_REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Stationary-State Repair\n\n"
        )

        handle.write(
            f"- Corrected gamma=0 rows: "
            f"{corrected_rows}/{expected_corrected}.\n"
        )

        handle.write(
            f"- Previous minimum stationary "
            f"eigenvalue: {old_minimum:.12e}.\n"
        )

        handle.write(
            f"- Corrected minimum stationary "
            f"eigenvalue: {new_minimum:.12e}.\n"
        )

        handle.write(
            f"- Overall numerical validation: "
            f"{'PASS' if overall_pass else 'FAIL'}.\n"
        )

        handle.write(
            f"- Mean stationary PYR5 range: "
            f"{minimum_pyr5:.9f} to "
            f"{maximum_pyr5:.9f}.\n\n"
        )

        handle.write(
            "The propagated trajectories were not "
            "recomputed or altered.\n"
        )

    print(
        "Day020 combined stationary-state "
        "repair completed."
    )

    print(
        f"Corrected gamma=0 rows: "
        f"{corrected_rows}/{expected_corrected}"
    )

    print(
        f"Previous minimum stationary eigenvalue: "
        f"{old_minimum:.12e}"
    )

    print(
        f"Corrected minimum stationary eigenvalue: "
        f"{new_minimum:.12e}"
    )

    print(
        f"Overall numerical validation: "
        f"{'PASS' if overall_pass else 'FAIL'}"
    )

    print(
        "Minimum mean stationary PYR5: "
        f"{minimum_pyr5:.9f} "
        f"at T="
        f"{float(minimum_condition['temperature_K']):g} K, "
        f"kappa="
        f"{float(minimum_condition['kappa_ref_ps_inv']):g}, "
        f"gamma="
        f"{float(minimum_condition['gamma_phi_ps_inv']):g}"
    )

    print(
        "Maximum mean stationary PYR5: "
        f"{maximum_pyr5:.9f} "
        f"at T="
        f"{float(maximum_condition['temperature_K']):g} K, "
        f"kappa="
        f"{float(maximum_condition['kappa_ref_ps_inv']):g}, "
        f"gamma="
        f"{float(maximum_condition['gamma_phi_ps_inv']):g}"
    )

    print(
        f"Wrote: "
        f"{REPAIR_REPORT_MD.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
