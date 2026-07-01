#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day020_confined_water_axial_radial_density"
)

OUTPUT_ROOT = (
    INPUT_ROOT
    / "regional_classification"
)

CANONICAL_NPZ = (
    INPUT_ROOT
    / "water_oxygen_axial_radial_density_canonical_v2.npz"
)

GEOMETRY_CSV = (
    INPUT_ROOT
    / "nanotube_geometry_summary.csv"
)

PYRENE_CSV = (
    INPUT_ROOT
    / "pyrene_geometry_summary.csv"
)

REGION_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "confined_water_region_summary.csv"
)

SENSITIVITY_CSV = (
    OUTPUT_ROOT
    / "confined_water_region_sensitivity.csv"
)

REGIONAL_NPZ = (
    OUTPUT_ROOT
    / "confined_water_regional_classification.npz"
)

FIGURE_STEM = (
    OUTPUT_ROOT
    / "figure_day020_confined_water_density_regions"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "CONFINED_WATER_REGIONAL_ANALYSIS_DAY020.md"
)

MANIFEST_MD = (
    OUTPUT_ROOT
    / "REGIONAL_ANALYSIS_MANIFEST_DAY020.md"
)

# Main operational definition.
INTERFACE_HALF_WIDTH_NM = 0.25
MOUTH_HALF_WIDTH_NM = 0.50

# Sensitivity grid.
INTERFACE_WIDTHS_NM = (
    0.20,
    0.25,
    0.30,
)

MOUTH_WIDTHS_NM = (
    0.40,
    0.50,
    0.60,
)

CONSERVATION_ABSOLUTE_TOLERANCE = 1.0e-8
DPI = 400

REGION_NAMES = (
    "interior_core",
    "interfacial_shell",
    "mouth_zones",
    "exterior",
)

REGION_CODES = {
    "interior_core": 1,
    "interfacial_shell": 2,
    "mouth_zones": 3,
    "exterior": 4,
}

REGION_DISPLAY_NAMES = {
    "interior_core": "Interior core",
    "interfacial_shell": "Interfacial shell",
    "mouth_zones": "Mouth zones",
    "exterior": "Exterior",
}


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def read_csv(
    path: Path,
) -> list[dict[str, str]]:
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


def load_geometry() -> dict[str, float]:
    rows = read_csv(
        GEOMETRY_CSV
    )

    if len(rows) != 1:
        raise RuntimeError(
            "Expected exactly one nanotube "
            "geometry row"
        )

    row = rows[0]

    required = (
        "axial_span_nm",
        "mean_wall_radius_nm",
        "p05_wall_radius_nm",
        "p95_wall_radius_nm",
    )

    missing = [
        key
        for key in required
        if key not in row
    ]

    if missing:
        raise RuntimeError(
            "Missing geometry fields: "
            + ", ".join(missing)
        )

    return {
        key: float(row[key])
        for key in required
    }


def load_pyrene_geometry(
) -> list[dict[str, float | str]]:
    rows = read_csv(
        PYRENE_CSV
    )

    if len(rows) != 4:
        raise RuntimeError(
            f"Expected four PYR residues, "
            f"found {len(rows)}"
        )

    result = []

    for row in rows:
        result.append(
            {
                "label": row[
                    "pyrene_label"
                ],
                "axial_position_nm": float(
                    row[
                        "axial_position_nm"
                    ]
                ),
                "radial_position_nm": float(
                    row[
                        "radial_position_nm"
                    ]
                ),
            }
        )

    return result


def load_density_data(
) -> dict[str, np.ndarray | float]:
    if not CANONICAL_NPZ.exists():
        raise RuntimeError(
            f"Missing canonical density product: "
            f"{CANONICAL_NPZ}"
        )

    with np.load(
        CANONICAL_NPZ,
        allow_pickle=False,
    ) as data:
        result = {
            key: data[key].copy()
            for key in data.files
        }

    required = (
        "axial_edges_nm",
        "radial_edges_nm",
        "axial_centers_nm",
        "radial_centers_nm",
        "cell_volumes_nm3",
        "density_nm3",
        "integrated_cell_counts",
        "integrated_water_count",
    )

    missing = [
        key
        for key in required
        if key not in result
    ]

    if missing:
        raise RuntimeError(
            "Canonical NPZ is missing: "
            + ", ".join(missing)
        )

    axial_edges = np.asarray(
        result["axial_edges_nm"],
        dtype=np.float64,
    )

    radial_edges = np.asarray(
        result["radial_edges_nm"],
        dtype=np.float64,
    )

    axial_centers = np.asarray(
        result["axial_centers_nm"],
        dtype=np.float64,
    )

    radial_centers = np.asarray(
        result["radial_centers_nm"],
        dtype=np.float64,
    )

    density = np.asarray(
        result["density_nm3"],
        dtype=np.float64,
    )

    volumes = np.asarray(
        result["cell_volumes_nm3"],
        dtype=np.float64,
    )

    counts = np.asarray(
        result["integrated_cell_counts"],
        dtype=np.float64,
    )

    expected_shape = (
        axial_centers.size,
        radial_centers.size,
    )

    for name, array in (
        ("density", density),
        ("cell volumes", volumes),
        ("integrated counts", counts),
    ):
        if array.shape != expected_shape:
            raise RuntimeError(
                f"Unexpected {name} shape: "
                f"{array.shape}; expected "
                f"{expected_shape}"
            )

    if axial_edges.size != (
        axial_centers.size + 1
    ):
        raise RuntimeError(
            "Invalid axial edge count"
        )

    if radial_edges.size != (
        radial_centers.size + 1
    ):
        raise RuntimeError(
            "Invalid radial edge count"
        )

    if not np.all(
        np.isfinite(density)
    ):
        raise RuntimeError(
            "Density contains non-finite values"
        )

    if np.min(density) < 0.0:
        raise RuntimeError(
            "Density contains negative values"
        )

    if np.min(volumes) <= 0.0:
        raise RuntimeError(
            "Cell volumes are not positive"
        )

    if np.min(counts) < 0.0:
        raise RuntimeError(
            "Integrated counts are negative"
        )

    return {
        "axial_edges_nm": axial_edges,
        "radial_edges_nm": radial_edges,
        "axial_centers_nm": axial_centers,
        "radial_centers_nm": radial_centers,
        "density_nm3": density,
        "cell_volumes_nm3": volumes,
        "integrated_cell_counts": counts,
        "integrated_water_count": float(
            np.asarray(
                result[
                    "integrated_water_count"
                ]
            )
        ),
    }


def classify_regions(
    axial_centers_nm: np.ndarray,
    radial_centers_nm: np.ndarray,
    tube_half_length_nm: float,
    wall_radius_nm: float,
    interface_half_width_nm: float,
    mouth_half_width_nm: float,
) -> dict[str, np.ndarray]:
    if interface_half_width_nm <= 0.0:
        raise RuntimeError(
            "Interface width must be positive"
        )

    if mouth_half_width_nm <= 0.0:
        raise RuntimeError(
            "Mouth width must be positive"
        )

    inner_radius_nm = (
        wall_radius_nm
        - interface_half_width_nm
    )

    outer_radius_nm = (
        wall_radius_nm
        + interface_half_width_nm
    )

    core_axial_limit_nm = (
        tube_half_length_nm
        - mouth_half_width_nm
    )

    outer_mouth_limit_nm = (
        tube_half_length_nm
        + mouth_half_width_nm
    )

    if inner_radius_nm <= 0.0:
        raise RuntimeError(
            "Interface width eliminates "
            "the interior radial region"
        )

    if core_axial_limit_nm <= 0.0:
        raise RuntimeError(
            "Mouth width eliminates "
            "the central axial region"
        )

    z_grid, r_grid = np.meshgrid(
        axial_centers_nm,
        radial_centers_nm,
        indexing="ij",
    )

    absolute_z = np.abs(
        z_grid
    )

    mouth_mask = (
        (absolute_z >= core_axial_limit_nm)
        & (
            absolute_z
            <= outer_mouth_limit_nm
        )
        & (
            r_grid < outer_radius_nm
        )
    )

    interfacial_mask = (
        (~mouth_mask)
        & (
            absolute_z
            < core_axial_limit_nm
        )
        & (
            r_grid
            >= inner_radius_nm
        )
        & (
            r_grid
            < outer_radius_nm
        )
    )

    interior_mask = (
        (~mouth_mask)
        & (~interfacial_mask)
        & (
            absolute_z
            < core_axial_limit_nm
        )
        & (
            r_grid
            < inner_radius_nm
        )
    )

    exterior_mask = ~(
        mouth_mask
        | interfacial_mask
        | interior_mask
    )

    masks = {
        "interior_core": interior_mask,
        "interfacial_shell": (
            interfacial_mask
        ),
        "mouth_zones": mouth_mask,
        "exterior": exterior_mask,
    }

    coverage = np.zeros_like(
        z_grid,
        dtype=np.int16,
    )

    for mask in masks.values():
        coverage += mask.astype(
            np.int16
        )

    if not np.all(
        coverage == 1
    ):
        unique, counts = np.unique(
            coverage,
            return_counts=True,
        )

        raise RuntimeError(
            "Regional masks do not form "
            "an exact partition: "
            f"{dict(zip(unique, counts))}"
        )

    return masks


def summarize_regions(
    masks: dict[str, np.ndarray],
    cell_volumes_nm3: np.ndarray,
    integrated_counts: np.ndarray,
    total_count: float,
    total_volume: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    summed_count = 0.0
    summed_volume = 0.0

    for region_name in REGION_NAMES:
        mask = masks[
            region_name
        ]

        region_count = float(
            np.sum(
                integrated_counts[
                    mask
                ]
            )
        )

        region_volume = float(
            np.sum(
                cell_volumes_nm3[
                    mask
                ]
            )
        )

        region_density = (
            region_count
            / region_volume
            if region_volume > 0.0
            else float("nan")
        )

        summed_count += region_count
        summed_volume += region_volume

        rows.append(
            {
                "region": region_name,
                "display_name": (
                    REGION_DISPLAY_NAMES[
                        region_name
                    ]
                ),
                "cell_count": int(
                    np.sum(mask)
                ),
                "volume_nm^3": (
                    region_volume
                ),
                "volume_fraction": (
                    region_volume
                    / total_volume
                ),
                "average_water_count": (
                    region_count
                ),
                "water_count_fraction": (
                    region_count
                    / total_count
                ),
                "volume_weighted_density_nm^-3": (
                    region_density
                ),
            }
        )

    count_error = abs(
        summed_count
        - total_count
    )

    volume_error = abs(
        summed_volume
        - total_volume
    )

    if count_error > (
        CONSERVATION_ABSOLUTE_TOLERANCE
    ):
        raise RuntimeError(
            "Regional water-count conservation "
            f"failed: error={count_error:.12e}"
        )

    if volume_error > (
        CONSERVATION_ABSOLUTE_TOLERANCE
    ):
        raise RuntimeError(
            "Regional volume conservation "
            f"failed: error={volume_error:.12e}"
        )

    return rows


def build_region_code_matrix(
    masks: dict[str, np.ndarray],
) -> np.ndarray:
    shape = next(
        iter(
            masks.values()
        )
    ).shape

    region_codes = np.zeros(
        shape,
        dtype=np.int8,
    )

    for region_name in REGION_NAMES:
        region_codes[
            masks[
                region_name
            ]
        ] = REGION_CODES[
            region_name
        ]

    if np.any(
        region_codes == 0
    ):
        raise RuntimeError(
            "Unclassified cells remain"
        )

    return region_codes


def axial_weighted_profile(
    density_nm3: np.ndarray,
    axial_edges_nm: np.ndarray,
    axial_mask: np.ndarray,
) -> np.ndarray:
    widths = np.diff(
        axial_edges_nm
    )

    selected_widths = widths[
        axial_mask
    ]

    if selected_widths.size == 0:
        raise RuntimeError(
            "No axial bins selected"
        )

    return np.sum(
        density_nm3[
            axial_mask,
            :
        ]
        * selected_widths[
            :,
            np.newaxis,
        ],
        axis=0,
    ) / np.sum(
        selected_widths
    )


def radial_weighted_profile(
    density_nm3: np.ndarray,
    radial_edges_nm: np.ndarray,
    radial_mask: np.ndarray,
) -> np.ndarray:
    annulus_areas = np.pi * (
        radial_edges_nm[1:] ** 2
        - radial_edges_nm[:-1] ** 2
    )

    selected_areas = annulus_areas[
        radial_mask
    ]

    if selected_areas.size == 0:
        raise RuntimeError(
            "No radial bins selected"
        )

    return np.sum(
        density_nm3[
            :,
            radial_mask,
        ]
        * selected_areas[
            np.newaxis,
            :
        ],
        axis=1,
    ) / np.sum(
        selected_areas
    )


def build_figure(
    density_data: dict[
        str,
        np.ndarray | float,
    ],
    geometry: dict[str, float],
    pyrene_rows: list[
        dict[str, float | str]
    ],
    main_summary: list[
        dict[str, object]
    ],
) -> None:
    axial_edges = np.asarray(
        density_data[
            "axial_edges_nm"
        ],
        dtype=np.float64,
    )

    radial_edges = np.asarray(
        density_data[
            "radial_edges_nm"
        ],
        dtype=np.float64,
    )

    axial_centers = np.asarray(
        density_data[
            "axial_centers_nm"
        ],
        dtype=np.float64,
    )

    radial_centers = np.asarray(
        density_data[
            "radial_centers_nm"
        ],
        dtype=np.float64,
    )

    density = np.asarray(
        density_data[
            "density_nm3"
        ],
        dtype=np.float64,
    )

    tube_half_length = (
        0.5
        * geometry[
            "axial_span_nm"
        ]
    )

    wall_radius = geometry[
        "mean_wall_radius_nm"
    ]

    inner_radius = (
        wall_radius
        - INTERFACE_HALF_WIDTH_NM
    )

    outer_radius = (
        wall_radius
        + INTERFACE_HALF_WIDTH_NM
    )

    core_axial_limit = (
        tube_half_length
        - MOUTH_HALF_WIDTH_NM
    )

    outer_mouth_limit = (
        tube_half_length
        + MOUTH_HALF_WIDTH_NM
    )

    central_axial_mask = (
        np.abs(
            axial_centers
        )
        < core_axial_limit
    )

    beyond_mouth_mask = (
        np.abs(
            axial_centers
        )
        > outer_mouth_limit
    )

    inner_radial_mask = (
        radial_centers
        < inner_radius
    )

    exterior_radial_mask = (
        radial_centers
        > outer_radius
    )

    central_radial_profile = (
        axial_weighted_profile(
            density,
            axial_edges,
            central_axial_mask,
        )
    )

    beyond_mouth_radial_profile = (
        axial_weighted_profile(
            density,
            axial_edges,
            beyond_mouth_mask,
        )
    )

    inner_axial_profile = (
        radial_weighted_profile(
            density,
            radial_edges,
            inner_radial_mask,
        )
    )

    exterior_axial_profile = (
        radial_weighted_profile(
            density,
            radial_edges,
            exterior_radial_mask,
        )
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.linewidth": 0.9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(13.2, 9.6),
    )

    axis_a = axes[0, 0]
    axis_b = axes[0, 1]
    axis_c = axes[1, 0]
    axis_d = axes[1, 1]

    mesh = axis_a.pcolormesh(
        axial_edges,
        radial_edges,
        density.T,
        shading="flat",
        cmap="viridis",
        vmin=0.0,
        vmax=float(
            np.max(density)
        ),
    )

    colorbar = figure.colorbar(
        mesh,
        ax=axis_a,
        pad=0.02,
    )

    colorbar.set_label(
        r"Water-O number density "
        r"(nm$^{-3}$)"
    )

    axis_a.axhline(
        wall_radius,
        linestyle="-",
        linewidth=1.4,
        color="white",
        label="Mean HBN wall",
    )

    axis_a.axhline(
        inner_radius,
        linestyle="--",
        linewidth=1.0,
        color="white",
    )

    axis_a.axhline(
        outer_radius,
        linestyle="--",
        linewidth=1.0,
        color="white",
    )

    for sign in (-1.0, 1.0):
        axis_a.axvline(
            sign * tube_half_length,
            linestyle="-",
            linewidth=1.3,
            color="white",
        )

        axis_a.axvline(
            sign * core_axial_limit,
            linestyle=":",
            linewidth=1.1,
            color="white",
        )

        axis_a.axvline(
            sign * outer_mouth_limit,
            linestyle=":",
            linewidth=1.1,
            color="white",
        )

    for pyrene in pyrene_rows:
        axis_a.scatter(
            float(
                pyrene[
                    "axial_position_nm"
                ]
            ),
            float(
                pyrene[
                    "radial_position_nm"
                ]
            ),
            marker="D",
            s=26,
            edgecolor="white",
            facecolor="black",
            linewidth=0.7,
            zorder=5,
        )

    axis_a.set_xlabel(
        "Axial coordinate z (nm)"
    )

    axis_a.set_ylabel(
        "Radial coordinate r (nm)"
    )

    axis_a.set_title(
        "Axial–radial water density"
    )

    axis_a.set_xlim(
        axial_edges[0],
        axial_edges[-1],
    )

    axis_a.set_ylim(
        radial_edges[0],
        radial_edges[-1],
    )

    axis_a.text(
        0.02,
        0.97,
        "(a)",
        transform=axis_a.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color="white",
    )

    axis_b.plot(
        radial_centers,
        central_radial_profile,
        label=(
            "Central tube span "
            r"($|z|<L/2-w_m$)"
        ),
    )

    axis_b.plot(
        radial_centers,
        beyond_mouth_radial_profile,
        label=(
            "Beyond both mouths "
            r"($|z|>L/2+w_m$)"
        ),
    )

    axis_b.axvline(
        wall_radius,
        linestyle="-",
        linewidth=1.2,
        color="black",
    )

    axis_b.axvspan(
        inner_radius,
        outer_radius,
        alpha=0.12,
    )

    axis_b.set_xlabel(
        "Radial coordinate r (nm)"
    )

    axis_b.set_ylabel(
        r"Axially averaged density "
        r"(nm$^{-3}$)"
    )

    axis_b.set_title(
        "Radial solvent organization"
    )

    axis_b.grid(
        True,
        alpha=0.22,
    )

    axis_b.legend(
        frameon=False,
    )

    axis_b.text(
        0.02,
        0.97,
        "(b)",
        transform=axis_b.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    axis_c.plot(
        axial_centers,
        inner_axial_profile,
        label=(
            r"Interior radial core "
            r"($r<R-w_i$)"
        ),
    )

    axis_c.plot(
        axial_centers,
        exterior_axial_profile,
        label=(
            r"Exterior radial region "
            r"($r>R+w_i$)"
        ),
    )

    for sign in (-1.0, 1.0):
        axis_c.axvline(
            sign * tube_half_length,
            linestyle="-",
            linewidth=1.1,
            color="black",
        )

        axis_c.axvspan(
            sign * core_axial_limit,
            sign * outer_mouth_limit,
            alpha=0.10,
        )

    axis_c.set_xlabel(
        "Axial coordinate z (nm)"
    )

    axis_c.set_ylabel(
        r"Radially averaged density "
        r"(nm$^{-3}$)"
    )

    axis_c.set_title(
        "Axial solvent organization"
    )

    axis_c.grid(
        True,
        alpha=0.22,
    )

    axis_c.legend(
        frameon=False,
    )

    axis_c.text(
        0.02,
        0.97,
        "(c)",
        transform=axis_c.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    display_names = [
        str(
            row[
                "display_name"
            ]
        )
        for row in main_summary
    ]

    counts = np.asarray(
        [
            float(
                row[
                    "average_water_count"
                ]
            )
            for row in main_summary
        ],
        dtype=np.float64,
    )

    fractions = np.asarray(
        [
            float(
                row[
                    "water_count_fraction"
                ]
            )
            for row in main_summary
        ],
        dtype=np.float64,
    )

    bars = axis_d.bar(
        np.arange(
            len(display_names)
        ),
        counts,
    )

    axis_d.set_xticks(
        np.arange(
            len(display_names)
        )
    )

    axis_d.set_xticklabels(
        display_names,
        rotation=18,
        ha="right",
    )

    axis_d.set_ylabel(
        "Average number of water molecules"
    )

    axis_d.set_title(
        "Regional population partition"
    )

    axis_d.grid(
        True,
        axis="y",
        alpha=0.22,
    )

    maximum_count = float(
        np.max(counts)
    )

    for bar, count, fraction in zip(
        bars,
        counts,
        fractions,
    ):
        axis_d.text(
            bar.get_x()
            + 0.5 * bar.get_width(),
            bar.get_height()
            + 0.018 * maximum_count,
            (
                f"{count:.1f}\n"
                f"({100.0 * fraction:.1f}%)"
            ),
            ha="center",
            va="bottom",
            fontsize=9,
        )

    axis_d.set_ylim(
        0.0,
        maximum_count * 1.18,
    )

    axis_d.text(
        0.02,
        0.97,
        "(d)",
        transform=axis_d.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    figure.suptitle(
        (
            "Confined-water organization around "
            "the frozen HBN–pyrene assembly"
        ),
        fontsize=14,
        y=0.985,
    )

    figure.tight_layout(
        rect=(
            0.0,
            0.0,
            1.0,
            0.965,
        )
    )

    figure.savefig(
        FIGURE_STEM.with_suffix(
            ".png"
        ),
        dpi=DPI,
        bbox_inches="tight",
    )

    figure.savefig(
        FIGURE_STEM.with_suffix(
            ".pdf"
        ),
        bbox_inches="tight",
    )

    plt.close(figure)


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    density_data = load_density_data()
    geometry = load_geometry()
    pyrene_rows = (
        load_pyrene_geometry()
    )

    axial_centers = np.asarray(
        density_data[
            "axial_centers_nm"
        ],
        dtype=np.float64,
    )

    radial_centers = np.asarray(
        density_data[
            "radial_centers_nm"
        ],
        dtype=np.float64,
    )

    cell_volumes = np.asarray(
        density_data[
            "cell_volumes_nm3"
        ],
        dtype=np.float64,
    )

    integrated_counts = np.asarray(
        density_data[
            "integrated_cell_counts"
        ],
        dtype=np.float64,
    )

    total_count = float(
        density_data[
            "integrated_water_count"
        ]
    )

    total_volume = float(
        np.sum(
            cell_volumes
        )
    )

    tube_half_length = (
        0.5
        * geometry[
            "axial_span_nm"
        ]
    )

    wall_radius = geometry[
        "mean_wall_radius_nm"
    ]

    main_masks = classify_regions(
        axial_centers,
        radial_centers,
        tube_half_length,
        wall_radius,
        INTERFACE_HALF_WIDTH_NM,
        MOUTH_HALF_WIDTH_NM,
    )

    main_summary = summarize_regions(
        main_masks,
        cell_volumes,
        integrated_counts,
        total_count,
        total_volume,
    )

    for row in main_summary:
        row[
            "interface_half_width_nm"
        ] = INTERFACE_HALF_WIDTH_NM

        row[
            "mouth_half_width_nm"
        ] = MOUTH_HALF_WIDTH_NM

        row[
            "tube_half_length_nm"
        ] = tube_half_length

        row[
            "wall_radius_nm"
        ] = wall_radius

    write_csv(
        REGION_SUMMARY_CSV,
        main_summary,
    )

    sensitivity_rows: list[
        dict[str, object]
    ] = []

    for interface_width in (
        INTERFACE_WIDTHS_NM
    ):
        for mouth_width in (
            MOUTH_WIDTHS_NM
        ):
            masks = classify_regions(
                axial_centers,
                radial_centers,
                tube_half_length,
                wall_radius,
                interface_width,
                mouth_width,
            )

            summary = summarize_regions(
                masks,
                cell_volumes,
                integrated_counts,
                total_count,
                total_volume,
            )

            for row in summary:
                sensitivity_rows.append(
                    {
                        "interface_half_width_nm": (
                            interface_width
                        ),
                        "mouth_half_width_nm": (
                            mouth_width
                        ),
                        **row,
                    }
                )

    write_csv(
        SENSITIVITY_CSV,
        sensitivity_rows,
    )

    region_codes = (
        build_region_code_matrix(
            main_masks
        )
    )

    np.savez_compressed(
        REGIONAL_NPZ,
        axial_centers_nm=np.asarray(
            density_data[
                "axial_centers_nm"
            ]
        ),
        radial_centers_nm=np.asarray(
            density_data[
                "radial_centers_nm"
            ]
        ),
        axial_edges_nm=np.asarray(
            density_data[
                "axial_edges_nm"
            ]
        ),
        radial_edges_nm=np.asarray(
            density_data[
                "radial_edges_nm"
            ]
        ),
        density_nm3=np.asarray(
            density_data[
                "density_nm3"
            ]
        ),
        cell_volumes_nm3=(
            cell_volumes
        ),
        integrated_cell_counts=(
            integrated_counts
        ),
        region_code=region_codes,
        region_names=np.asarray(
            REGION_NAMES
        ),
        region_code_values=np.asarray(
            [
                REGION_CODES[name]
                for name in REGION_NAMES
            ],
            dtype=np.int8,
        ),
        tube_half_length_nm=np.asarray(
            tube_half_length
        ),
        wall_radius_nm=np.asarray(
            wall_radius
        ),
        interface_half_width_nm=np.asarray(
            INTERFACE_HALF_WIDTH_NM
        ),
        mouth_half_width_nm=np.asarray(
            MOUTH_HALF_WIDTH_NM
        ),
        total_integrated_water_count=np.asarray(
            total_count
        ),
    )

    build_figure(
        density_data,
        geometry,
        pyrene_rows,
        main_summary,
    )

    count_values = {
        str(row["region"]): float(
            row[
                "average_water_count"
            ]
        )
        for row in main_summary
    }

    fraction_values = {
        str(row["region"]): float(
            row[
                "water_count_fraction"
            ]
        )
        for row in main_summary
    }

    reconstructed_count = sum(
        count_values.values()
    )

    conservation_error = abs(
        reconstructed_count
        - total_count
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Confined-Water "
            "Regional Analysis\n\n"
        )

        handle.write(
            "## Scope\n\n"
        )

        handle.write(
            "The time-averaged axial–radial "
            "water-oxygen density was integrated "
            "using the exact cylindrical volume "
            "of every grid cell. The analyzed "
            "domain contains an average of "
            f"{total_count:.6f} water molecules.\n\n"
        )

        handle.write(
            "## Operational region definition\n\n"
        )

        handle.write(
            f"- Nanotube half-length: "
            f"{tube_half_length:.6f} nm.\n"
        )

        handle.write(
            f"- Mean HBN wall radius: "
            f"{wall_radius:.6f} nm.\n"
        )

        handle.write(
            f"- Interfacial half-width: "
            f"{INTERFACE_HALF_WIDTH_NM:.3f} nm.\n"
        )

        handle.write(
            f"- Mouth half-width: "
            f"{MOUTH_HALF_WIDTH_NM:.3f} nm.\n\n"
        )

        handle.write(
            "The four regions are mutually "
            "exclusive and cover every grid cell. "
            "Classification is based on grid-cell "
            "centers, while populations use exact "
            "cylindrical cell volumes.\n\n"
        )

        handle.write(
            "## Main regional populations\n\n"
        )

        for region_name in REGION_NAMES:
            handle.write(
                f"- "
                f"{REGION_DISPLAY_NAMES[region_name]}: "
                f"{count_values[region_name]:.6f} "
                "molecules "
                f"({100.0 * fraction_values[region_name]:.3f}%).\n"
            )

        handle.write(
            "\n## Conservation\n\n"
        )

        handle.write(
            f"- Reconstructed population: "
            f"{reconstructed_count:.9f}.\n"
        )

        handle.write(
            f"- Canonical population: "
            f"{total_count:.9f}.\n"
        )

        handle.write(
            f"- Absolute conservation error: "
            f"{conservation_error:.3e}.\n\n"
        )

        handle.write(
            "## Sensitivity analysis\n\n"
        )

        handle.write(
            "The regional partition was repeated "
            "for interfacial half-widths of "
            "0.20, 0.25, and 0.30 nm and mouth "
            "half-widths of 0.40, 0.50, and "
            "0.60 nm. These operational choices "
            "change the numerical allocation among "
            "neighboring regions but do not alter "
            "the underlying density field.\n\n"
        )

        handle.write(
            "## Interpretation boundary\n\n"
        )

        handle.write(
            "These are time-averaged occupancies "
            "around a frozen solute. They describe "
            "solvent organization but not coupled "
            "water–solute conformational dynamics, "
            "water residence times, or scaffold "
            "thermal stability.\n"
        )

    with MANIFEST_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Confined-Water "
            "Regional Analysis Manifest\n\n"
        )

        for path in (
            REGION_SUMMARY_CSV,
            SENSITIVITY_CSV,
            REGIONAL_NPZ,
            FIGURE_STEM.with_suffix(
                ".png"
            ),
            FIGURE_STEM.with_suffix(
                ".pdf"
            ),
            REPORT_MD,
        ):
            handle.write(
                f"- `{relative(path)}`\n"
            )

    required_outputs = (
        REGION_SUMMARY_CSV,
        SENSITIVITY_CSV,
        REGIONAL_NPZ,
        FIGURE_STEM.with_suffix(
            ".png"
        ),
        FIGURE_STEM.with_suffix(
            ".pdf"
        ),
        REPORT_MD,
        MANIFEST_MD,
    )

    missing = [
        path
        for path in required_outputs
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty outputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    log(
        "Day020 confined-water regional "
        "analysis completed."
    )

    log(
        f"Tube half-length: "
        f"{tube_half_length:.6f} nm"
    )

    log(
        f"Mean wall radius: "
        f"{wall_radius:.6f} nm"
    )

    log(
        f"Interface half-width: "
        f"{INTERFACE_HALF_WIDTH_NM:.3f} nm"
    )

    log(
        f"Mouth half-width: "
        f"{MOUTH_HALF_WIDTH_NM:.3f} nm"
    )

    log(
        f"Canonical cylinder count: "
        f"{total_count:.6f}"
    )

    for region_name in REGION_NAMES:
        log(
            f"{REGION_DISPLAY_NAMES[region_name]}: "
            f"{count_values[region_name]:.6f} "
            f"({100.0 * fraction_values[region_name]:.3f}%)"
        )

    log(
        f"Regional conservation error: "
        f"{conservation_error:.3e}"
    )

    log(
        f"Sensitivity conditions: "
        f"{len(INTERFACE_WIDTHS_NM) * len(MOUTH_WIDTHS_NM)}"
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
