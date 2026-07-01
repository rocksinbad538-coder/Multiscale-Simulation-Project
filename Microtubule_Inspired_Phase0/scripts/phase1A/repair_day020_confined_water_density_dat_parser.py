#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day020_confined_water_axial_radial_density"
)

DENSITY_DAT = (
    INPUT_ROOT
    / "water_oxygen_axial_radial_density.dat"
)

PARAMETERS_CSV = (
    INPUT_ROOT
    / "densmap_parameters.csv"
)

GEOMETRY_CSV = (
    INPUT_ROOT
    / "nanotube_geometry_summary.csv"
)

CANONICAL_NPZ = (
    INPUT_ROOT
    / "water_oxygen_axial_radial_density_canonical.npz"
)

VALIDATION_CSV = (
    INPUT_ROOT
    / "density_dat_parser_validation.csv"
)

REPORT_MD = (
    INPUT_ROOT
    / "DENSITY_DAT_PARSER_REPAIR_DAY020.md"
)

NEGATIVE_DENSITY_TOLERANCE = 1.0e-10


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def read_key_value_csv(
    path: Path,
) -> dict[str, str]:
    if not path.exists():
        return {}

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        return {}

    if (
        "parameter" in rows[0]
        and "value" in rows[0]
    ):
        return {
            row["parameter"]: row["value"]
            for row in rows
        }

    return rows[0]


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


def parse_density_dat() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    if (
        not DENSITY_DAT.exists()
        or DENSITY_DAT.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty DAT file: "
            f"{DENSITY_DAT}"
        )

    raw = np.loadtxt(
        DENSITY_DAT,
        dtype=np.float64,
    )

    if raw.ndim != 2:
        raise RuntimeError(
            f"Expected a 2D DAT matrix, "
            f"found ndim={raw.ndim}"
        )

    if (
        raw.shape[0] < 2
        or raw.shape[1] < 2
    ):
        raise RuntimeError(
            f"DAT matrix is too small: "
            f"{raw.shape}"
        )

    radial_centers_nm = raw[
        0,
        1:,
    ].copy()

    axial_centers_nm = raw[
        1:,
        0,
    ].copy()

    density_nm3 = raw[
        1:,
        1:,
    ].copy()

    if density_nm3.shape != (
        axial_centers_nm.size,
        radial_centers_nm.size,
    ):
        raise RuntimeError(
            "Coordinate and density dimensions "
            "are inconsistent"
        )

    return (
        raw,
        axial_centers_nm,
        radial_centers_nm,
        density_nm3,
    )


def validate_coordinates(
    axial_centers_nm: np.ndarray,
    radial_centers_nm: np.ndarray,
) -> tuple[float, float]:
    if not np.all(
        np.isfinite(
            axial_centers_nm
        )
    ):
        raise RuntimeError(
            "Axial coordinates contain "
            "non-finite values"
        )

    if not np.all(
        np.isfinite(
            radial_centers_nm
        )
    ):
        raise RuntimeError(
            "Radial coordinates contain "
            "non-finite values"
        )

    axial_differences = np.diff(
        axial_centers_nm
    )

    radial_differences = np.diff(
        radial_centers_nm
    )

    if not np.all(
        axial_differences > 0.0
    ):
        raise RuntimeError(
            "Axial coordinates are not "
            "strictly increasing"
        )

    if not np.all(
        radial_differences > 0.0
    ):
        raise RuntimeError(
            "Radial coordinates are not "
            "strictly increasing"
        )

    if np.min(
        radial_centers_nm
    ) < -1.0e-12:
        raise RuntimeError(
            "Radial coordinates contain "
            "negative values"
        )

    axial_spacing_nm = float(
        np.median(
            axial_differences
        )
    )

    radial_spacing_nm = float(
        np.median(
            radial_differences
        )
    )

    return (
        axial_spacing_nm,
        radial_spacing_nm,
    )


def validate_density(
    density_nm3: np.ndarray,
) -> tuple[int, float]:
    if not np.all(
        np.isfinite(
            density_nm3
        )
    ):
        raise RuntimeError(
            "Density matrix contains "
            "non-finite values"
        )

    negative_mask = (
        density_nm3
        < -NEGATIVE_DENSITY_TOLERANCE
    )

    significant_negative_count = int(
        np.sum(
            negative_mask
        )
    )

    if significant_negative_count > 0:
        raise RuntimeError(
            "Density matrix contains "
            f"{significant_negative_count} "
            "significantly negative values"
        )

    tiny_negative_count = int(
        np.sum(
            (
                density_nm3 < 0.0
            )
            & (
                density_nm3
                >= -NEGATIVE_DENSITY_TOLERANCE
            )
        )
    )

    if tiny_negative_count > 0:
        density_nm3[
            density_nm3 < 0.0
        ] = 0.0

    return (
        tiny_negative_count,
        float(
            np.min(
                density_nm3
            )
        ),
    )


def coordinate_edges(
    centers: np.ndarray,
) -> np.ndarray:
    if centers.size < 2:
        raise RuntimeError(
            "At least two coordinate centers "
            "are required"
        )

    internal_edges = 0.5 * (
        centers[:-1]
        + centers[1:]
    )

    first_edge = (
        centers[0]
        - 0.5
        * (
            centers[1]
            - centers[0]
        )
    )

    last_edge = (
        centers[-1]
        + 0.5
        * (
            centers[-1]
            - centers[-2]
        )
    )

    return np.concatenate(
        (
            [first_edge],
            internal_edges,
            [last_edge],
        )
    )


def main() -> None:
    (
        raw,
        axial_centers_nm,
        radial_centers_nm,
        density_nm3,
    ) = parse_density_dat()

    (
        axial_spacing_nm,
        radial_spacing_nm,
    ) = validate_coordinates(
        axial_centers_nm,
        radial_centers_nm,
    )

    (
        tiny_negative_density_count,
        minimum_density_nm3,
    ) = validate_density(
        density_nm3
    )

    axial_edges_nm = coordinate_edges(
        axial_centers_nm
    )

    radial_edges_nm = coordinate_edges(
        radial_centers_nm
    )

    negative_raw_entries = int(
        np.sum(
            raw < 0.0
        )
    )

    negative_axial_coordinates = int(
        np.sum(
            axial_centers_nm < 0.0
        )
    )

    negative_radial_coordinates = int(
        np.sum(
            radial_centers_nm < 0.0
        )
    )

    negative_density_entries = int(
        np.sum(
            density_nm3 < 0.0
        )
    )

    if negative_raw_entries != (
        negative_axial_coordinates
    ):
        raise RuntimeError(
            "Negative raw entries are not "
            "fully explained by axial coordinates"
        )

    if negative_radial_coordinates != 0:
        raise RuntimeError(
            "Unexpected negative radial coordinates"
        )

    if negative_density_entries != 0:
        raise RuntimeError(
            "Negative density values remain "
            "after validation"
        )

    parameters = read_key_value_csv(
        PARAMETERS_CSV
    )

    geometry = read_key_value_csv(
        GEOMETRY_CSV
    )

    np.savez_compressed(
        CANONICAL_NPZ,
        axial_centers_nm=(
            axial_centers_nm
        ),
        radial_centers_nm=(
            radial_centers_nm
        ),
        axial_edges_nm=(
            axial_edges_nm
        ),
        radial_edges_nm=(
            radial_edges_nm
        ),
        density_nm3=density_nm3,
        density_unit=np.asarray(
            "nm^-3"
        ),
        source_dat=np.asarray(
            relative(
                DENSITY_DAT
            )
        ),
    )

    validation_rows = [
        {
            "source_DAT": relative(
                DENSITY_DAT
            ),
            "raw_rows": raw.shape[0],
            "raw_columns": raw.shape[1],
            "axial_bin_count": (
                axial_centers_nm.size
            ),
            "radial_bin_count": (
                radial_centers_nm.size
            ),
            "density_rows": (
                density_nm3.shape[0]
            ),
            "density_columns": (
                density_nm3.shape[1]
            ),
            "density_cells": (
                density_nm3.size
            ),
            "axial_min_center_nm": float(
                axial_centers_nm[0]
            ),
            "axial_max_center_nm": float(
                axial_centers_nm[-1]
            ),
            "radial_min_center_nm": float(
                radial_centers_nm[0]
            ),
            "radial_max_center_nm": float(
                radial_centers_nm[-1]
            ),
            "median_axial_spacing_nm": (
                axial_spacing_nm
            ),
            "median_radial_spacing_nm": (
                radial_spacing_nm
            ),
            "raw_negative_entries": (
                negative_raw_entries
            ),
            "negative_axial_coordinates": (
                negative_axial_coordinates
            ),
            "negative_radial_coordinates": (
                negative_radial_coordinates
            ),
            "negative_density_entries": (
                negative_density_entries
            ),
            "tiny_negative_density_values_clipped": (
                tiny_negative_density_count
            ),
            "minimum_density_nm^-3": (
                minimum_density_nm3
            ),
            "maximum_density_nm^-3": float(
                np.max(
                    density_nm3
                )
            ),
            "unweighted_grid_mean_density_nm^-3": float(
                np.mean(
                    density_nm3
                )
            ),
            "finite_coordinates_pass": True,
            "monotonic_coordinates_pass": True,
            "finite_density_pass": True,
            "nonnegative_density_pass": True,
            "coordinate_header_interpretation_pass": True,
            "overall_validation_pass": True,
        }
    ]

    write_csv(
        VALIDATION_CSV,
        validation_rows,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Density DAT Parser Repair\n\n"
        )

        handle.write(
            "## Identified DAT structure\n\n"
        )

        handle.write(
            "The GROMACS textual density output "
            "contains coordinate headers:\n\n"
        )

        handle.write(
            "- Row 0, columns 1 onward: "
            "radial-bin centers.\n"
        )

        handle.write(
            "- Column 0, rows 1 onward: "
            "axial-bin centers.\n"
        )

        handle.write(
            "- Rows 1 onward, columns 1 onward: "
            "physical water-oxygen number density.\n\n"
        )

        handle.write(
            "The negative entries detected by the "
            "original validator were axial coordinates "
            "at negative z, not negative densities.\n\n"
        )

        handle.write(
            "## Canonical grid\n\n"
        )

        handle.write(
            f"- Raw DAT dimensions: "
            f"{raw.shape[0]} × "
            f"{raw.shape[1]}.\n"
        )

        handle.write(
            f"- Density dimensions: "
            f"{density_nm3.shape[0]} × "
            f"{density_nm3.shape[1]}.\n"
        )

        handle.write(
            f"- Axial centers: "
            f"{axial_centers_nm[0]:.6f} to "
            f"{axial_centers_nm[-1]:.6f} nm.\n"
        )

        handle.write(
            f"- Radial centers: "
            f"{radial_centers_nm[0]:.6f} to "
            f"{radial_centers_nm[-1]:.6f} nm.\n"
        )

        handle.write(
            f"- Median axial spacing: "
            f"{axial_spacing_nm:.8f} nm.\n"
        )

        handle.write(
            f"- Median radial spacing: "
            f"{radial_spacing_nm:.8f} nm.\n"
        )

        handle.write(
            f"- Minimum density: "
            f"{np.min(density_nm3):.6f} nm^-3.\n"
        )

        handle.write(
            f"- Maximum density: "
            f"{np.max(density_nm3):.6f} nm^-3.\n\n"
        )

        handle.write(
            "## Cross-checks against analysis inputs\n\n"
        )

        handle.write(
            f"- Requested bin width: "
            f"{parameters.get('bin_width_nm', 'not available')} nm.\n"
        )

        handle.write(
            f"- Requested axial half-range: "
            f"{parameters.get('axial_half_range_nm', 'not available')} nm.\n"
        )

        handle.write(
            f"- Requested radial maximum: "
            f"{parameters.get('radial_maximum_nm', 'not available')} nm.\n"
        )

        handle.write(
            f"- HBN axial span: "
            f"{geometry.get('axial_span_nm', 'not available')} nm.\n"
        )

        handle.write(
            f"- Mean HBN wall radius: "
            f"{geometry.get('mean_wall_radius_nm', 'not available')} nm.\n\n"
        )

        handle.write(
            "## Canonical product\n\n"
        )

        handle.write(
            f"- `{relative(CANONICAL_NPZ)}`\n\n"
        )

        handle.write(
            "The NPZ file stores axial and radial "
            "centers, bin edges, and the validated "
            "160 × 54 density matrix. It is the "
            "canonical input for subsequent plotting "
            "and regional solvent classification.\n\n"
        )

        handle.write(
            "## Status\n\n"
        )

        handle.write(
            "The existing density calculation is "
            "numerically valid. No GROMACS rerun and "
            "no XPM conversion are required.\n"
        )

    log(
        "Day020 density DAT parser repair completed."
    )

    log(
        f"Raw DAT shape: "
        f"{raw.shape[0]} x {raw.shape[1]}"
    )

    log(
        f"Canonical density shape: "
        f"{density_nm3.shape[0]} x "
        f"{density_nm3.shape[1]}"
    )

    log(
        f"Axial range: "
        f"{axial_centers_nm[0]:.6f} to "
        f"{axial_centers_nm[-1]:.6f} nm"
    )

    log(
        f"Radial range: "
        f"{radial_centers_nm[0]:.6f} to "
        f"{radial_centers_nm[-1]:.6f} nm"
    )

    log(
        "Raw negative entries: "
        f"{negative_raw_entries}"
    )

    log(
        "Negative axial coordinates: "
        f"{negative_axial_coordinates}"
    )

    log(
        "Negative density entries: "
        f"{negative_density_entries}"
    )

    log(
        "Density range: "
        f"{np.min(density_nm3):.6f} to "
        f"{np.max(density_nm3):.6f} nm^-3"
    )

    log(
        "Overall validation: PASS"
    )

    log(
        f"Canonical NPZ: "
        f"{relative(CANONICAL_NPZ)}"
    )

    log(
        f"Wrote: "
        f"{relative(INPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
