#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ACCEPTED_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day020_confined_water_axial_radial_density"
)

XTC_PATH = (
    ACCEPTED_ROOT
    / "nvt_100ps_frozenSolute.xtc"
)

TPR_PATH = (
    ACCEPTED_ROOT
    / "nvt_100ps_frozenSolute.tpr"
)

INDEX_PATH = (
    OUTPUT_ROOT
    / "confined_water_analysis_groups.ndx"
)

PARAMETERS_CSV = (
    OUTPUT_ROOT
    / "densmap_parameters.csv"
)

DENSITY_DAT = (
    OUTPUT_ROOT
    / "water_oxygen_axial_radial_density.dat"
)

DENSITY_XPM = (
    OUTPUT_ROOT
    / "water_oxygen_axial_radial_density.xpm"
)

XPM_LOG = (
    OUTPUT_ROOT
    / "gromacs_densmap_xpm_only.log"
)

VALIDATION_CSV = (
    OUTPUT_ROOT
    / "density_output_validation.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "CONFINED_WATER_DENSITY_FINALIZATION_DAY020.md"
)


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def find_gromacs() -> str:
    candidates = []

    configured = os.environ.get("GMX_BIN")

    if configured:
        candidates.append(configured)

    candidates.extend(
        [
            "/usr/local/gromacs/bin/gmx",
            "gmx",
            "gmx_mpi",
        ]
    )

    for candidate in candidates:
        resolved = shutil.which(candidate)

        if resolved is not None:
            return resolved

    raise RuntimeError(
        "Could not locate GROMACS. "
        "Set GMX_BIN explicitly."
    )


def read_parameters() -> dict[str, str]:
    if not PARAMETERS_CSV.exists():
        raise RuntimeError(
            f"Missing parameter file: "
            f"{PARAMETERS_CSV}"
        )

    with PARAMETERS_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    parameters = {
        row["parameter"]: row["value"]
        for row in rows
    }

    required = (
        "begin_time_ps",
        "end_time_ps",
        "bin_width_nm",
        "axial_half_range_nm",
        "radial_maximum_nm",
        "density_unit",
    )

    missing = [
        key
        for key in required
        if key not in parameters
    ]

    if missing:
        raise RuntimeError(
            "Missing densmap parameters: "
            + ", ".join(missing)
        )

    return parameters


def validate_dat() -> tuple[
    np.ndarray,
    int,
    int,
]:
    if (
        not DENSITY_DAT.exists()
        or DENSITY_DAT.stat().st_size == 0
    ):
        raise RuntimeError(
            "The existing density DAT file "
            "is missing or empty."
        )

    matrix = np.loadtxt(
        DENSITY_DAT,
        dtype=np.float64,
    )

    if matrix.ndim == 1:
        matrix = matrix[np.newaxis, :]

    if matrix.ndim != 2:
        raise RuntimeError(
            f"Unexpected DAT dimensionality: "
            f"{matrix.ndim}"
        )

    if matrix.size == 0:
        raise RuntimeError(
            "Density DAT matrix is empty."
        )

    if not np.all(
        np.isfinite(matrix)
    ):
        raise RuntimeError(
            "Density DAT contains "
            "non-finite values."
        )

    if np.min(matrix) < -1.0e-12:
        raise RuntimeError(
            "Density DAT contains "
            "negative density values."
        )

    return (
        matrix,
        int(matrix.shape[0]),
        int(matrix.shape[1]),
    )


def run_xpm_only(
    gromacs: str,
    parameters: dict[str, str],
) -> None:
    command = [
        gromacs,
        "densmap",
        "-f",
        str(XTC_PATH),
        "-s",
        str(TPR_PATH),
        "-n",
        str(INDEX_PATH),
        "-o",
        str(DENSITY_XPM),
        "-b",
        parameters[
            "begin_time_ps"
        ],
        "-e",
        parameters[
            "end_time_ps"
        ],
        "-bin",
        parameters[
            "bin_width_nm"
        ],
        "-amax",
        parameters[
            "axial_half_range_nm"
        ],
        "-rmax",
        parameters[
            "radial_maximum_nm"
        ],
        "-unit",
        "nm-3",
    ]

    log("Generating auxiliary XPM output.")
    log(" ".join(command))
    log()

    with XPM_LOG.open(
        "w",
        encoding="utf-8",
    ) as log_handle:
        process = subprocess.Popen(
            command,
            cwd=PROJECT_ROOT,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        if (
            process.stdin is None
            or process.stdout is None
        ):
            raise RuntimeError(
                "Could not open densmap streams."
            )

        process.stdin.write(
            "AxisMinus\n"
            "AxisPlus\n"
            "Water_O\n"
        )

        process.stdin.close()

        for line in process.stdout:
            print(
                line,
                end="",
                flush=True,
            )

            log_handle.write(line)
            log_handle.flush()

        return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            f"XPM-only densmap execution "
            f"failed with code {return_code}"
        )

    if (
        not DENSITY_XPM.exists()
        or DENSITY_XPM.stat().st_size == 0
    ):
        raise RuntimeError(
            "The XPM-only execution finished, "
            "but the XPM file is missing or empty."
        )


def write_validation(
    matrix: np.ndarray,
    n_rows: int,
    n_columns: int,
) -> None:
    rows = [
        {
            "DAT_path": relative(
                DENSITY_DAT
            ),
            "DAT_size_bytes": (
                DENSITY_DAT.stat().st_size
            ),
            "XPM_path": relative(
                DENSITY_XPM
            ),
            "XPM_size_bytes": (
                DENSITY_XPM.stat().st_size
            ),
            "DAT_rows": n_rows,
            "DAT_columns": n_columns,
            "DAT_cells": int(
                matrix.size
            ),
            "minimum_density_nm^-3": float(
                np.min(matrix)
            ),
            "maximum_density_nm^-3": float(
                np.max(matrix)
            ),
            "mean_grid_density_nm^-3": float(
                np.mean(matrix)
            ),
            "finite_values_pass": bool(
                np.all(
                    np.isfinite(matrix)
                )
            ),
            "nonnegative_density_pass": bool(
                np.min(matrix)
                >= -1.0e-12
            ),
            "DAT_validation_pass": True,
            "XPM_validation_pass": True,
            "overall_validation_pass": True,
        }
    ]

    with VALIDATION_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                rows[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(rows)


def write_report(
    matrix: np.ndarray,
    n_rows: int,
    n_columns: int,
    xpm_was_generated: bool,
) -> None:
    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Confined-Water "
            "Density Output Finalization\n\n"
        )

        handle.write(
            "## Correction\n\n"
        )

        handle.write(
            "The original GROMACS execution "
            "successfully generated the quantitative "
            "DAT density matrix. The script then "
            "incorrectly required an XPM product from "
            "the same execution. GROMACS treats `-od` "
            "and `-o` as alternative output modes.\n\n"
        )

        handle.write(
            "The existing DAT result was preserved. "
            "A separate XPM-only execution was used "
            "to generate the visualization matrix.\n\n"
        )

        handle.write(
            "## Validated products\n\n"
        )

        handle.write(
            f"- DAT matrix: "
            f"`{relative(DENSITY_DAT)}`.\n"
        )

        handle.write(
            f"- XPM matrix: "
            f"`{relative(DENSITY_XPM)}`.\n"
        )

        handle.write(
            f"- DAT dimensions: "
            f"{n_rows} × {n_columns}.\n"
        )

        handle.write(
            f"- DAT cells: "
            f"{matrix.size}.\n"
        )

        handle.write(
            f"- Minimum grid density: "
            f"{float(np.min(matrix)):.9f} nm^-3.\n"
        )

        handle.write(
            f"- Maximum grid density: "
            f"{float(np.max(matrix)):.9f} nm^-3.\n"
        )

        handle.write(
            f"- Mean grid density: "
            f"{float(np.mean(matrix)):.9f} nm^-3.\n"
        )

        handle.write(
            f"- XPM generated in this step: "
            f"{xpm_was_generated}.\n\n"
        )

        handle.write(
            "## Status\n\n"
        )

        handle.write(
            "The axial–radial confined-water "
            "density calculation is complete and "
            "numerically valid. The DAT matrix is "
            "the canonical quantitative product; "
            "the XPM file is an auxiliary "
            "visualization product.\n"
        )


def main() -> None:
    required_inputs = (
        XTC_PATH,
        TPR_PATH,
        INDEX_PATH,
        PARAMETERS_CSV,
        DENSITY_DAT,
    )

    for path in required_inputs:
        if not path.exists():
            raise RuntimeError(
                f"Missing required input: {path}"
            )

    parameters = read_parameters()

    (
        matrix,
        n_rows,
        n_columns,
    ) = validate_dat()

    xpm_was_generated = False

    if (
        not DENSITY_XPM.exists()
        or DENSITY_XPM.stat().st_size == 0
    ):
        run_xpm_only(
            find_gromacs(),
            parameters,
        )

        xpm_was_generated = True

    write_validation(
        matrix,
        n_rows,
        n_columns,
    )

    write_report(
        matrix,
        n_rows,
        n_columns,
        xpm_was_generated,
    )

    log()
    log(
        "Day020 confined-water density "
        "finalization completed."
    )

    log(
        f"DAT dimensions: "
        f"{n_rows} x {n_columns}"
    )

    log(
        "Density range: "
        f"{float(np.min(matrix)):.6f} to "
        f"{float(np.max(matrix)):.6f} nm^-3"
    )

    log(
        f"XPM generated now: "
        f"{xpm_was_generated}"
    )

    log(
        "Overall validation: PASS"
    )

    log(
        f"Wrote: "
        f"{relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
