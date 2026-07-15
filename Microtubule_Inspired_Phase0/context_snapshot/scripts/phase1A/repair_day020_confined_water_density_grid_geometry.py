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

CANONICAL_V1 = (
    INPUT_ROOT
    / "water_oxygen_axial_radial_density_canonical.npz"
)

CANONICAL_V2 = (
    INPUT_ROOT
    / "water_oxygen_axial_radial_density_canonical_v2.npz"
)

VALIDATION_CSV = (
    INPUT_ROOT
    / "density_grid_geometry_validation.csv"
)

REPORT_MD = (
    INPUT_ROOT
    / "DENSITY_GRID_GEOMETRY_REPAIR_DAY020.md"
)

CURRENT_POINTER = (
    INPUT_ROOT
    / "CANONICAL_DENSITY_PRODUCT.txt"
)

COORDINATE_TOLERANCE_NM = 5.0e-5
SPACING_RELATIVE_TOLERANCE = 1.0e-3
NEGATIVE_DENSITY_TOLERANCE = 1.0e-12


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def read_parameter_csv() -> dict[str, str]:
    if not PARAMETERS_CSV.exists():
        raise RuntimeError(
            f"Missing parameters CSV: "
            f"{PARAMETERS_CSV}"
        )

    with PARAMETERS_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    parameters = {
        row["parameter"]: row["value"]
        for row in rows
    }

    required = (
        "axial_half_range_nm",
        "radial_maximum_nm",
        "bin_width_nm",
        "begin_time_ps",
        "end_time_ps",
    )

    missing = [
        name
        for name in required
        if name not in parameters
    ]

    if missing:
        raise RuntimeError(
            "Missing parameters: "
            + ", ".join(missing)
        )

    return parameters


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


def spacing_statistics(
    lower_edges: np.ndarray,
) -> tuple[float, float, float]:
    differences = np.diff(
        lower_edges
    )

    if not np.all(
        differences > 0.0
    ):
        raise RuntimeError(
            "Grid coordinates are not "
            "strictly increasing"
        )

    median_spacing = float(
        np.median(
            differences
        )
    )

    minimum_spacing = float(
        np.min(
            differences
        )
    )

    maximum_spacing = float(
        np.max(
            differences
        )
    )

    relative_deviation = float(
        np.max(
            np.abs(
                differences
                - median_spacing
            )
        )
        / median_spacing
    )

    if relative_deviation > (
        SPACING_RELATIVE_TOLERANCE
    ):
        raise RuntimeError(
            "Grid spacing is insufficiently "
            "uniform: relative deviation "
            f"{relative_deviation:.6e}"
        )

    return (
        median_spacing,
        minimum_spacing,
        maximum_spacing,
    )


def main() -> None:
    if (
        not DENSITY_DAT.exists()
        or DENSITY_DAT.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty density DAT: "
            f"{DENSITY_DAT}"
        )

    parameters = read_parameter_csv()

    axial_half_range_nm = float(
        parameters[
            "axial_half_range_nm"
        ]
    )

    radial_maximum_nm = float(
        parameters[
            "radial_maximum_nm"
        ]
    )

    requested_bin_width_nm = float(
        parameters[
            "bin_width_nm"
        ]
    )

    raw = np.loadtxt(
        DENSITY_DAT,
        dtype=np.float64,
    )

    if raw.ndim != 2:
        raise RuntimeError(
            f"Expected 2-D DAT matrix, "
            f"found ndim={raw.ndim}"
        )

    if raw.shape != (161, 55):
        raise RuntimeError(
            f"Unexpected DAT shape: "
            f"{raw.shape}"
        )

    axial_lower_edges_nm = raw[
        1:,
        0,
    ].copy()

    radial_lower_edges_nm = raw[
        0,
        1:,
    ].copy()

    density_nm3 = raw[
        1:,
        1:,
    ].copy()

    if density_nm3.shape != (
        axial_lower_edges_nm.size,
        radial_lower_edges_nm.size,
    ):
        raise RuntimeError(
            "Coordinate and density dimensions "
            "are inconsistent"
        )

    if not np.all(
        np.isfinite(
            density_nm3
        )
    ):
        raise RuntimeError(
            "Density matrix contains "
            "non-finite values"
        )

    significant_negative_density = int(
        np.sum(
            density_nm3
            < -NEGATIVE_DENSITY_TOLERANCE
        )
    )

    if significant_negative_density:
        raise RuntimeError(
            "Physical density matrix contains "
            f"{significant_negative_density} "
            "negative entries"
        )

    tiny_negative_mask = (
        density_nm3 < 0.0
    )

    tiny_negative_count = int(
        np.sum(
            tiny_negative_mask
        )
    )

    if tiny_negative_count:
        density_nm3[
            tiny_negative_mask
        ] = 0.0

    (
        axial_nominal_spacing_nm,
        axial_minimum_spacing_nm,
        axial_maximum_spacing_nm,
    ) = spacing_statistics(
        axial_lower_edges_nm
    )

    (
        radial_nominal_spacing_nm,
        radial_minimum_spacing_nm,
        radial_maximum_spacing_nm,
    ) = spacing_statistics(
        radial_lower_edges_nm
    )

    axial_start_error_nm = abs(
        axial_lower_edges_nm[0]
        + axial_half_range_nm
    )

    radial_start_error_nm = abs(
        radial_lower_edges_nm[0]
    )

    predicted_axial_upper_nm = (
        axial_lower_edges_nm[-1]
        + axial_nominal_spacing_nm
    )

    predicted_radial_upper_nm = (
        radial_lower_edges_nm[-1]
        + radial_nominal_spacing_nm
    )

    axial_upper_error_nm = abs(
        predicted_axial_upper_nm
        - axial_half_range_nm
    )

    radial_upper_error_nm = abs(
        predicted_radial_upper_nm
        - radial_maximum_nm
    )

    coordinate_errors = {
        "axial_start_error_nm": (
            axial_start_error_nm
        ),
        "radial_start_error_nm": (
            radial_start_error_nm
        ),
        "axial_upper_error_nm": (
            axial_upper_error_nm
        ),
        "radial_upper_error_nm": (
            radial_upper_error_nm
        ),
    }

    failed_errors = {
        name: error
        for name, error
        in coordinate_errors.items()
        if error > COORDINATE_TOLERANCE_NM
    }

    if failed_errors:
        raise RuntimeError(
            "Coordinate-boundary interpretation "
            f"failed: {failed_errors}"
        )

    axial_edges_nm = np.concatenate(
        (
            axial_lower_edges_nm,
            np.asarray(
                [
                    axial_half_range_nm
                ],
                dtype=np.float64,
            ),
        )
    )

    radial_edges_nm = np.concatenate(
        (
            radial_lower_edges_nm,
            np.asarray(
                [
                    radial_maximum_nm
                ],
                dtype=np.float64,
            ),
        )
    )

    if not np.all(
        np.diff(
            axial_edges_nm
        ) > 0.0
    ):
        raise RuntimeError(
            "Final axial edges are invalid"
        )

    if not np.all(
        np.diff(
            radial_edges_nm
        ) > 0.0
    ):
        raise RuntimeError(
            "Final radial edges are invalid"
        )

    if radial_edges_nm[0] != 0.0:
        raise RuntimeError(
            "Physical radial grid does not "
            "start at zero"
        )

    axial_centers_nm = 0.5 * (
        axial_edges_nm[:-1]
        + axial_edges_nm[1:]
    )

    radial_centers_nm = 0.5 * (
        radial_edges_nm[:-1]
        + radial_edges_nm[1:]
    )

    axial_widths_nm = np.diff(
        axial_edges_nm
    )

    radial_annulus_areas_nm2 = (
        np.pi
        * (
            radial_edges_nm[1:] ** 2
            - radial_edges_nm[:-1] ** 2
        )
    )

    cell_volumes_nm3 = (
        axial_widths_nm[:, np.newaxis]
        * radial_annulus_areas_nm2[
            np.newaxis,
            :
        ]
    )

    if cell_volumes_nm3.shape != (
        density_nm3.shape
    ):
        raise RuntimeError(
            "Cell-volume and density grids "
            "have different shapes"
        )

    if not np.all(
        cell_volumes_nm3 > 0.0
    ):
        raise RuntimeError(
            "Non-positive cylindrical "
            "cell volume detected"
        )

    integrated_cell_counts = (
        density_nm3
        * cell_volumes_nm3
    )

    integrated_water_count = float(
        np.sum(
            integrated_cell_counts
        )
    )

    analyzed_cylinder_volume_nm3 = float(
        np.sum(
            cell_volumes_nm3
        )
    )

    volume_weighted_density_nm3 = (
        integrated_water_count
        / analyzed_cylinder_volume_nm3
    )

    direct_cylinder_volume_nm3 = (
        np.pi
        * radial_maximum_nm**2
        * (
            2.0
            * axial_half_range_nm
        )
    )

    cylinder_volume_error_nm3 = abs(
        analyzed_cylinder_volume_nm3
        - direct_cylinder_volume_nm3
    )

    cylinder_volume_relative_error = (
        cylinder_volume_error_nm3
        / direct_cylinder_volume_nm3
    )

    if cylinder_volume_relative_error > (
        1.0e-12
    ):
        raise RuntimeError(
            "Cylindrical cell volumes do not "
            "recover the analytical domain volume"
        )

    np.savez_compressed(
        CANONICAL_V2,
        axial_lower_edges_nm=(
            axial_lower_edges_nm
        ),
        radial_lower_edges_nm=(
            radial_lower_edges_nm
        ),
        axial_edges_nm=(
            axial_edges_nm
        ),
        radial_edges_nm=(
            radial_edges_nm
        ),
        axial_centers_nm=(
            axial_centers_nm
        ),
        radial_centers_nm=(
            radial_centers_nm
        ),
        axial_widths_nm=(
            axial_widths_nm
        ),
        radial_annulus_areas_nm2=(
            radial_annulus_areas_nm2
        ),
        cell_volumes_nm3=(
            cell_volumes_nm3
        ),
        density_nm3=(
            density_nm3
        ),
        integrated_cell_counts=(
            integrated_cell_counts
        ),
        integrated_water_count=np.asarray(
            integrated_water_count
        ),
        analyzed_cylinder_volume_nm3=np.asarray(
            analyzed_cylinder_volume_nm3
        ),
        volume_weighted_density_nm3=np.asarray(
            volume_weighted_density_nm3
        ),
        axial_half_range_nm=np.asarray(
            axial_half_range_nm
        ),
        radial_maximum_nm=np.asarray(
            radial_maximum_nm
        ),
        requested_bin_width_nm=np.asarray(
            requested_bin_width_nm
        ),
        density_unit=np.asarray(
            "nm^-3"
        ),
        source_dat=np.asarray(
            relative(
                DENSITY_DAT
            )
        ),
        coordinate_semantics=np.asarray(
            "DAT row/column coordinates are "
            "lower cell edges"
        ),
    )

    validation_rows = [
        {
            "source_DAT": relative(
                DENSITY_DAT
            ),
            "superseded_canonical_v1": (
                relative(CANONICAL_V1)
                if CANONICAL_V1.exists()
                else "not_found"
            ),
            "canonical_v2": relative(
                CANONICAL_V2
            ),
            "raw_rows": raw.shape[0],
            "raw_columns": raw.shape[1],
            "density_rows": (
                density_nm3.shape[0]
            ),
            "density_columns": (
                density_nm3.shape[1]
            ),
            "axial_edge_count": (
                axial_edges_nm.size
            ),
            "radial_edge_count": (
                radial_edges_nm.size
            ),
            "axial_center_count": (
                axial_centers_nm.size
            ),
            "radial_center_count": (
                radial_centers_nm.size
            ),
            "axial_min_edge_nm": (
                axial_edges_nm[0]
            ),
            "axial_max_edge_nm": (
                axial_edges_nm[-1]
            ),
            "radial_min_edge_nm": (
                radial_edges_nm[0]
            ),
            "radial_max_edge_nm": (
                radial_edges_nm[-1]
            ),
            "axial_nominal_spacing_nm": (
                axial_nominal_spacing_nm
            ),
            "radial_nominal_spacing_nm": (
                radial_nominal_spacing_nm
            ),
            "requested_bin_width_nm": (
                requested_bin_width_nm
            ),
            "minimum_density_nm^-3": float(
                np.min(
                    density_nm3
                )
            ),
            "maximum_density_nm^-3": float(
                np.max(
                    density_nm3
                )
            ),
            "tiny_negative_values_clipped": (
                tiny_negative_count
            ),
            "analyzed_cylinder_volume_nm^3": (
                analyzed_cylinder_volume_nm3
            ),
            "integrated_average_water_count": (
                integrated_water_count
            ),
            "volume_weighted_density_nm^-3": (
                volume_weighted_density_nm3
            ),
            "cylinder_volume_relative_error": (
                cylinder_volume_relative_error
            ),
            "lower_edge_interpretation_pass": True,
            "physical_radial_origin_pass": True,
            "density_nonnegative_pass": True,
            "cylindrical_volume_pass": True,
            "overall_validation_pass": True,
        }
    ]

    write_csv(
        VALIDATION_CSV,
        validation_rows,
    )

    CURRENT_POINTER.write_text(
        (
            "Canonical quantitative density product:\n"
            f"{relative(CANONICAL_V2)}\n\n"
            "The previous canonical NPZ interpreted "
            "DAT coordinates as cell centers and is "
            "retained only for traceability.\n"
        ),
        encoding="utf-8",
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Density Grid Geometry Repair\n\n"
        )

        handle.write(
            "## Correct coordinate semantics\n\n"
        )

        handle.write(
            "The first DAT column and row contain "
            "lower cell edges, not cell centers. "
            "This follows from their exact agreement "
            "with the requested domain boundaries:\n\n"
        )

        handle.write(
            f"- Axial start: "
            f"{axial_lower_edges_nm[0]:.6f} nm.\n"
        )

        handle.write(
            f"- Requested axial domain: "
            f"[-{axial_half_range_nm:.6f}, "
            f"+{axial_half_range_nm:.6f}] nm.\n"
        )

        handle.write(
            f"- Radial start: "
            f"{radial_lower_edges_nm[0]:.6f} nm.\n"
        )

        handle.write(
            f"- Requested radial domain: "
            f"[0, {radial_maximum_nm:.6f}] nm.\n\n"
        )

        handle.write(
            "The physical bin centers were therefore "
            "reconstructed from consecutive edges. "
            "The radial grid now starts at exactly "
            "zero and contains no artificial negative "
            "radial boundary.\n\n"
        )

        handle.write(
            "## Corrected grid\n\n"
        )

        handle.write(
            f"- Density cells: "
            f"{density_nm3.shape[0]} × "
            f"{density_nm3.shape[1]}.\n"
        )

        handle.write(
            f"- Axial edges: "
            f"{axial_edges_nm[0]:.6f} to "
            f"{axial_edges_nm[-1]:.6f} nm.\n"
        )

        handle.write(
            f"- Radial edges: "
            f"{radial_edges_nm[0]:.6f} to "
            f"{radial_edges_nm[-1]:.6f} nm.\n"
        )

        handle.write(
            f"- Median axial spacing: "
            f"{axial_nominal_spacing_nm:.8f} nm.\n"
        )

        handle.write(
            f"- Median radial spacing: "
            f"{radial_nominal_spacing_nm:.8f} nm.\n\n"
        )

        handle.write(
            "## Cylindrical integration\n\n"
        )

        handle.write(
            "Each grid cell was assigned the exact "
            "cylindrical volume\n\n"
        )

        handle.write(
            "\\[\n"
            "\\Delta V_{ij}="
            "\\pi\\left(r_{j+1}^{2}-r_j^{2}\\right)"
            "\\left(z_{i+1}-z_i\\right).\n"
            "\\]\n\n"
        )

        handle.write(
            f"- Analysis-domain volume: "
            f"{analyzed_cylinder_volume_nm3:.9f} nm³.\n"
        )

        handle.write(
            f"- Integrated average water count in "
            f"the analyzed cylinder: "
            f"{integrated_water_count:.9f}.\n"
        )

        handle.write(
            f"- Volume-weighted mean density: "
            f"{volume_weighted_density_nm3:.9f} nm⁻³.\n\n"
        )

        handle.write(
            "## Product status\n\n"
        )

        handle.write(
            f"- Current canonical product: "
            f"`{relative(CANONICAL_V2)}`.\n"
        )

        handle.write(
            f"- Previous product retained for "
            f"traceability: "
            f"`{relative(CANONICAL_V1)}`.\n\n"
        )

        handle.write(
            "The V2 product must be used for all "
            "figures and regional cylindrical "
            "integrations.\n"
        )

    log(
        "Day020 density-grid geometry repair completed."
    )

    log(
        f"Raw DAT shape: "
        f"{raw.shape[0]} x {raw.shape[1]}"
    )

    log(
        f"Density grid: "
        f"{density_nm3.shape[0]} x "
        f"{density_nm3.shape[1]}"
    )

    log(
        "Coordinate interpretation: "
        "lower cell edges"
    )

    log(
        "Axial edges: "
        f"{axial_edges_nm[0]:.6f} to "
        f"{axial_edges_nm[-1]:.6f} nm"
    )

    log(
        "Radial edges: "
        f"{radial_edges_nm[0]:.6f} to "
        f"{radial_edges_nm[-1]:.6f} nm"
    )

    log(
        "Axial nominal spacing: "
        f"{axial_nominal_spacing_nm:.8f} nm"
    )

    log(
        "Radial nominal spacing: "
        f"{radial_nominal_spacing_nm:.8f} nm"
    )

    log(
        "Density range: "
        f"{np.min(density_nm3):.6f} to "
        f"{np.max(density_nm3):.6f} nm^-3"
    )

    log(
        "Analysis cylinder volume: "
        f"{analyzed_cylinder_volume_nm3:.6f} nm^3"
    )

    log(
        "Integrated average water count: "
        f"{integrated_water_count:.6f}"
    )

    log(
        "Volume-weighted mean density: "
        f"{volume_weighted_density_nm3:.6f} nm^-3"
    )

    log(
        "Overall validation: PASS"
    )

    log(
        f"Canonical V2: "
        f"{relative(CANONICAL_V2)}"
    )


if __name__ == "__main__":
    main()
