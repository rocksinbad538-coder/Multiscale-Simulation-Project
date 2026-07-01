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
    / "profile_guided_classification"
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

HBN_SEGMENTS_CSV = (
    INPUT_ROOT
    / "hbn_architecture_audit/"
    "hbn_axial_segments.csv"
)

HBN_SUMMARY_CSV = (
    INPUT_ROOT
    / "hbn_architecture_audit/"
    "hbn_architecture_summary.csv"
)

OPERATIONAL_SUMMARY_CSV = (
    INPUT_ROOT
    / "regional_classification/"
    "confined_water_region_summary.csv"
)

BOUNDARIES_CSV = (
    OUTPUT_ROOT
    / "profile_guided_boundaries.csv"
)

REGION_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "profile_guided_region_summary.csv"
)

SENSITIVITY_CSV = (
    OUTPUT_ROOT
    / "profile_guided_region_sensitivity.csv"
)

COMPARISON_CSV = (
    OUTPUT_ROOT
    / "operational_vs_profile_guided_regions.csv"
)

RADIAL_PROFILE_CSV = (
    OUTPUT_ROOT
    / "profile_guided_radial_profile.csv"
)

AXIAL_PROFILE_CSV = (
    OUTPUT_ROOT
    / "profile_guided_axial_profile.csv"
)

REGIONAL_NPZ = (
    OUTPUT_ROOT
    / "profile_guided_regional_classification.npz"
)

FIGURE_STEM = (
    OUTPUT_ROOT
    / "figure_day020_confined_water_profile_guided_regions"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "CONFINED_WATER_PROFILE_GUIDED_ANALYSIS_DAY020.md"
)

MANIFEST_MD = (
    OUTPUT_ROOT
    / "PROFILE_GUIDED_ANALYSIS_MANIFEST_DAY020.md"
)

# Profile construction.
END_EXCLUSION_NM = 0.50
RADIAL_SMOOTHING_SIGMA_NM = 0.10
AXIAL_SMOOTHING_SIGMA_NM = 0.15

# Radial depletion-shell definition.
MAIN_RADIAL_RECOVERY_FRACTION = 0.50
RADIAL_RECOVERY_FRACTIONS = (
    0.40,
    0.50,
    0.60,
)

# Axial mouth-transition definition.
MAIN_AXIAL_INNER_FRACTION = 0.80
MAIN_AXIAL_OUTER_FRACTION = 0.20

AXIAL_INNER_FRACTIONS = (
    0.75,
    0.80,
    0.85,
)

AXIAL_SEARCH_HALF_WIDTH_NM = 1.00
OUTSIDE_PLATEAU_OFFSET_NM = 0.35
CENTRAL_PLATEAU_FRACTION = 0.35

MINIMUM_RADIAL_CONTRAST_NM3 = 3.0
MINIMUM_AXIAL_CONTRAST_NM3 = 3.0

CONSERVATION_TOLERANCE = 1.0e-8
DPI = 400

REGION_NAMES = (
    "lumen_core",
    "interfacial_shell",
    "mouth_transitions",
    "exterior",
)

REGION_LABELS = {
    "lumen_core": "Lumen core",
    "interfacial_shell": "Interfacial shell",
    "mouth_transitions": "Mouth transitions",
    "exterior": "Exterior",
}

REGION_CODES = {
    "lumen_core": 1,
    "interfacial_shell": 2,
    "mouth_transitions": 3,
    "exterior": 4,
}


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


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

    fieldnames: list[str] = []
    seen_fields: set[str] = set()

    for row in rows:
        for fieldname in row.keys():
            if fieldname not in seen_fields:
                fieldnames.append(
                    fieldname
                )
                seen_fields.add(
                    fieldname
                )

    normalized_rows = [
        {
            fieldname: row.get(
                fieldname,
                "",
            )
            for fieldname in fieldnames
        }
        for row in rows
    ]

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(
            normalized_rows
        )


def gaussian_smooth(
    values: np.ndarray,
    spacing: float,
    sigma: float,
) -> np.ndarray:
    if sigma <= 0.0:
        return values.copy()

    radius = max(
        2,
        int(
            np.ceil(
                4.0 * sigma / spacing
            )
        ),
    )

    offsets = np.arange(
        -radius,
        radius + 1,
        dtype=np.float64,
    ) * spacing

    kernel = np.exp(
        -0.5
        * (
            offsets / sigma
        ) ** 2
    )

    kernel /= np.sum(kernel)

    padded = np.pad(
        values,
        radius,
        mode="edge",
    )

    smoothed = np.convolve(
        padded,
        kernel,
        mode="same",
    )

    return smoothed[
        radius:-radius
    ]


def interpolate_crossing(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    target: float,
) -> float:
    if np.isclose(y1, y0):
        return 0.5 * (
            x0 + x1
        )

    fraction = (
        target - y0
    ) / (
        y1 - y0
    )

    fraction = float(
        np.clip(
            fraction,
            0.0,
            1.0,
        )
    )

    return x0 + fraction * (
        x1 - x0
    )


def find_crossings(
    coordinates: np.ndarray,
    values: np.ndarray,
    target: float,
    mask: np.ndarray,
) -> list[float]:
    indices = np.where(mask)[0]

    if indices.size < 2:
        return []

    crossings: list[float] = []

    for left, right in zip(
        indices[:-1],
        indices[1:],
    ):
        if right != left + 1:
            continue

        y0 = values[left] - target
        y1 = values[right] - target

        if y0 == 0.0:
            crossings.append(
                float(
                    coordinates[left]
                )
            )
            continue

        if y1 == 0.0:
            crossings.append(
                float(
                    coordinates[right]
                )
            )
            continue

        if y0 * y1 < 0.0:
            crossings.append(
                interpolate_crossing(
                    float(
                        coordinates[left]
                    ),
                    float(
                        values[left]
                    ),
                    float(
                        coordinates[right]
                    ),
                    float(
                        values[right]
                    ),
                    target,
                )
            )

    return crossings


def nearest_crossing(
    crossings: list[float],
    reference: float,
    label: str,
) -> float:
    if not crossings:
        raise RuntimeError(
            f"No crossing found for {label}"
        )

    return min(
        crossings,
        key=lambda value: abs(
            value - reference
        ),
    )


def load_density() -> dict[str, np.ndarray | float]:
    if not CANONICAL_NPZ.exists():
        raise RuntimeError(
            f"Missing canonical density: "
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
        "density_nm3",
        "cell_volumes_nm3",
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
            "Canonical density product is missing: "
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
        ("volumes", volumes),
        ("counts", counts),
    ):
        if array.shape != expected_shape:
            raise RuntimeError(
                f"Unexpected {name} shape: "
                f"{array.shape}"
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


def load_hbn_geometry() -> dict[str, float]:
    segment_rows = read_csv(
        HBN_SEGMENTS_CSV
    )

    if len(segment_rows) != 1:
        raise RuntimeError(
            "Profile-guided analysis requires "
            "one continuous HBN segment"
        )

    summary_rows = read_csv(
        HBN_SUMMARY_CSV
    )

    if len(summary_rows) != 1:
        raise RuntimeError(
            "Expected one HBN architecture summary row"
        )

    segment = segment_rows[0]
    summary = summary_rows[0]

    if int(
        summary[
            "detected_segment_count"
        ]
    ) != 1:
        raise RuntimeError(
            "HBN architecture is not continuous"
        )

    return {
        "lower_boundary_nm": float(
            segment[
                "lower_boundary_nm"
            ]
        ),
        "upper_boundary_nm": float(
            segment[
                "upper_boundary_nm"
            ]
        ),
        "mean_radius_nm": float(
            segment[
                "mean_radius_nm"
            ]
        ),
        "p05_radius_nm": float(
            segment[
                "p05_radius_nm"
            ]
        ),
        "p95_radius_nm": float(
            segment[
                "p95_radius_nm"
            ]
        ),
    }


def load_pyrenes(
) -> list[dict[str, float | str]]:
    rows = read_csv(
        PYRENE_CSV
    )

    if len(rows) != 4:
        raise RuntimeError(
            f"Expected four PYR rows, "
            f"found {len(rows)}"
        )

    return [
        {
            "label": row[
                "pyrene_label"
            ],
            "z_nm": float(
                row[
                    "axial_position_nm"
                ]
            ),
            "r_nm": float(
                row[
                    "radial_position_nm"
                ]
            ),
        }
        for row in rows
    ]


def axial_average_radial_profile(
    density: np.ndarray,
    axial_edges: np.ndarray,
    axial_mask: np.ndarray,
) -> np.ndarray:
    widths = np.diff(
        axial_edges
    )

    selected_widths = widths[
        axial_mask
    ]

    if selected_widths.size == 0:
        raise RuntimeError(
            "No axial bins selected "
            "for radial profile"
        )

    return (
        np.sum(
            density[
                axial_mask,
                :
            ]
            * selected_widths[
                :,
                np.newaxis,
            ],
            axis=0,
        )
        / np.sum(
            selected_widths
        )
    )


def radial_average_axial_profile(
    density: np.ndarray,
    radial_edges: np.ndarray,
    radial_mask: np.ndarray,
) -> np.ndarray:
    annulus_areas = np.pi * (
        radial_edges[1:] ** 2
        - radial_edges[:-1] ** 2
    )

    selected_areas = annulus_areas[
        radial_mask
    ]

    if selected_areas.size == 0:
        raise RuntimeError(
            "No radial bins selected "
            "for axial profile"
        )

    return (
        np.sum(
            density[
                :,
                radial_mask,
            ]
            * selected_areas[
                np.newaxis,
                :
            ],
            axis=1,
        )
        / np.sum(
            selected_areas
        )
    )


def derive_radial_boundaries(
    radial_centers: np.ndarray,
    radial_profile_smooth: np.ndarray,
    wall_radius: float,
    recovery_fraction: float,
) -> dict[str, float]:
    search_mask = (
        radial_centers
        >= wall_radius - 0.40
    ) & (
        radial_centers
        <= wall_radius + 0.40
    )

    search_indices = np.where(
        search_mask
    )[0]

    if search_indices.size < 3:
        raise RuntimeError(
            "Insufficient radial bins "
            "around HBN wall"
        )

    minimum_index = int(
        search_indices[
            np.argmin(
                radial_profile_smooth[
                    search_indices
                ]
            )
        ]
    )

    minimum_radius = float(
        radial_centers[
            minimum_index
        ]
    )

    minimum_density = float(
        radial_profile_smooth[
            minimum_index
        ]
    )

    inner_plateau_mask = (
        radial_centers >= 0.20
    ) & (
        radial_centers
        <= wall_radius - 0.45
    )

    outer_plateau_mask = (
        radial_centers
        >= wall_radius + 0.55
    ) & (
        radial_centers
        <= min(
            radial_centers[-1] - 0.15,
            wall_radius + 1.20,
        )
    )

    if np.sum(
        inner_plateau_mask
    ) < 3:
        raise RuntimeError(
            "Insufficient inner plateau bins"
        )

    if np.sum(
        outer_plateau_mask
    ) < 3:
        raise RuntimeError(
            "Insufficient outer plateau bins"
        )

    inner_plateau_density = float(
        np.median(
            radial_profile_smooth[
                inner_plateau_mask
            ]
        )
    )

    outer_plateau_density = float(
        np.median(
            radial_profile_smooth[
                outer_plateau_mask
            ]
        )
    )

    inner_contrast = (
        inner_plateau_density
        - minimum_density
    )

    outer_contrast = (
        outer_plateau_density
        - minimum_density
    )

    if inner_contrast < (
        MINIMUM_RADIAL_CONTRAST_NM3
    ):
        raise RuntimeError(
            "Inner radial depletion contrast "
            f"is too small: {inner_contrast:.6f}"
        )

    if outer_contrast < (
        MINIMUM_RADIAL_CONTRAST_NM3
    ):
        raise RuntimeError(
            "Outer radial depletion contrast "
            f"is too small: {outer_contrast:.6f}"
        )

    inner_target = (
        minimum_density
        + recovery_fraction
        * inner_contrast
    )

    outer_target = (
        minimum_density
        + recovery_fraction
        * outer_contrast
    )

    inner_search_mask = (
        radial_centers
        <= minimum_radius
    )

    outer_search_mask = (
        radial_centers
        >= minimum_radius
    )

    inner_crossings = find_crossings(
        radial_centers,
        radial_profile_smooth,
        inner_target,
        inner_search_mask,
    )

    outer_crossings = find_crossings(
        radial_centers,
        radial_profile_smooth,
        outer_target,
        outer_search_mask,
    )

    inner_boundary = nearest_crossing(
        inner_crossings,
        minimum_radius,
        "inner radial boundary",
    )

    outer_boundary = nearest_crossing(
        outer_crossings,
        minimum_radius,
        "outer radial boundary",
    )

    if not (
        inner_boundary
        < minimum_radius
        < outer_boundary
    ):
        raise RuntimeError(
            "Derived radial boundaries "
            "do not bracket the minimum"
        )

    return {
        "recovery_fraction": (
            recovery_fraction
        ),
        "minimum_radius_nm": (
            minimum_radius
        ),
        "minimum_density_nm^-3": (
            minimum_density
        ),
        "inner_plateau_density_nm^-3": (
            inner_plateau_density
        ),
        "outer_plateau_density_nm^-3": (
            outer_plateau_density
        ),
        "inner_target_density_nm^-3": (
            inner_target
        ),
        "outer_target_density_nm^-3": (
            outer_target
        ),
        "inner_boundary_nm": (
            inner_boundary
        ),
        "outer_boundary_nm": (
            outer_boundary
        ),
        "interfacial_width_nm": (
            outer_boundary
            - inner_boundary
        ),
    }


def derive_axial_boundaries(
    axial_centers: np.ndarray,
    axial_profile_smooth: np.ndarray,
    segment_lower: float,
    segment_upper: float,
    inner_fraction: float,
    outer_fraction: float,
) -> dict[str, float]:
    segment_center = 0.5 * (
        segment_lower
        + segment_upper
    )

    segment_half_span = 0.5 * (
        segment_upper
        - segment_lower
    )

    central_mask = (
        np.abs(
            axial_centers
            - segment_center
        )
        <= (
            CENTRAL_PLATEAU_FRACTION
            * segment_half_span
        )
    )

    left_outside_mask = (
        axial_centers
        <= (
            segment_lower
            - OUTSIDE_PLATEAU_OFFSET_NM
        )
    )

    right_outside_mask = (
        axial_centers
        >= (
            segment_upper
            + OUTSIDE_PLATEAU_OFFSET_NM
        )
    )

    if np.sum(
        central_mask
    ) < 5:
        raise RuntimeError(
            "Insufficient central axial bins"
        )

    if np.sum(
        left_outside_mask
    ) < 3:
        raise RuntimeError(
            "Insufficient left outside bins"
        )

    if np.sum(
        right_outside_mask
    ) < 3:
        raise RuntimeError(
            "Insufficient right outside bins"
        )

    central_density = float(
        np.median(
            axial_profile_smooth[
                central_mask
            ]
        )
    )

    left_outside_density = float(
        np.median(
            axial_profile_smooth[
                left_outside_mask
            ]
        )
    )

    right_outside_density = float(
        np.median(
            axial_profile_smooth[
                right_outside_mask
            ]
        )
    )

    left_contrast = (
        central_density
        - left_outside_density
    )

    right_contrast = (
        central_density
        - right_outside_density
    )

    if left_contrast < (
        MINIMUM_AXIAL_CONTRAST_NM3
    ):
        raise RuntimeError(
            "Left mouth axial contrast "
            f"is too small: {left_contrast:.6f}"
        )

    if right_contrast < (
        MINIMUM_AXIAL_CONTRAST_NM3
    ):
        raise RuntimeError(
            "Right mouth axial contrast "
            f"is too small: {right_contrast:.6f}"
        )

    left_inner_target = (
        left_outside_density
        + inner_fraction
        * left_contrast
    )

    left_outer_target = (
        left_outside_density
        + outer_fraction
        * left_contrast
    )

    right_inner_target = (
        right_outside_density
        + inner_fraction
        * right_contrast
    )

    right_outer_target = (
        right_outside_density
        + outer_fraction
        * right_contrast
    )

    # Search each transition over its complete axial
    # half-domain. The previous fixed ±1 nm window could
    # exclude high-recovery crossings when the solvent
    # transition extends farther into the lumen.
    left_search_mask = (
        axial_centers
        <= segment_center
    )

    right_search_mask = (
        axial_centers
        >= segment_center
    )

    left_inner = nearest_crossing(
        find_crossings(
            axial_centers,
            axial_profile_smooth,
            left_inner_target,
            left_search_mask,
        ),
        segment_lower,
        "left inner mouth boundary",
    )

    left_outer = nearest_crossing(
        find_crossings(
            axial_centers,
            axial_profile_smooth,
            left_outer_target,
            left_search_mask,
        ),
        segment_lower,
        "left outer mouth boundary",
    )

    right_inner = nearest_crossing(
        find_crossings(
            axial_centers,
            axial_profile_smooth,
            right_inner_target,
            right_search_mask,
        ),
        segment_upper,
        "right inner mouth boundary",
    )

    right_outer = nearest_crossing(
        find_crossings(
            axial_centers,
            axial_profile_smooth,
            right_outer_target,
            right_search_mask,
        ),
        segment_upper,
        "right outer mouth boundary",
    )

    # Left side increases from outside to interior.
    if left_outer > left_inner:
        left_outer, left_inner = (
            left_inner,
            left_outer,
        )

    # Right side decreases from interior to outside.
    if right_inner > right_outer:
        right_inner, right_outer = (
            right_outer,
            right_inner,
        )

    if not (
        left_outer
        < left_inner
        < right_inner
        < right_outer
    ):
        raise RuntimeError(
            "Derived axial mouth boundaries "
            "are not correctly ordered"
        )

    return {
        "inner_fraction": (
            inner_fraction
        ),
        "outer_fraction": (
            outer_fraction
        ),
        "central_density_nm^-3": (
            central_density
        ),
        "left_outside_density_nm^-3": (
            left_outside_density
        ),
        "right_outside_density_nm^-3": (
            right_outside_density
        ),
        "left_inner_target_nm^-3": (
            left_inner_target
        ),
        "left_outer_target_nm^-3": (
            left_outer_target
        ),
        "right_inner_target_nm^-3": (
            right_inner_target
        ),
        "right_outer_target_nm^-3": (
            right_outer_target
        ),
        "left_outer_boundary_nm": (
            left_outer
        ),
        "left_inner_boundary_nm": (
            left_inner
        ),
        "right_inner_boundary_nm": (
            right_inner
        ),
        "right_outer_boundary_nm": (
            right_outer
        ),
        "left_transition_width_nm": (
            left_inner
            - left_outer
        ),
        "right_transition_width_nm": (
            right_outer
            - right_inner
        ),
    }


def classify_regions(
    axial_centers: np.ndarray,
    radial_centers: np.ndarray,
    radial_boundaries: dict[str, float],
    axial_boundaries: dict[str, float],
) -> dict[str, np.ndarray]:
    z_grid, r_grid = np.meshgrid(
        axial_centers,
        radial_centers,
        indexing="ij",
    )

    inner_r = radial_boundaries[
        "inner_boundary_nm"
    ]

    outer_r = radial_boundaries[
        "outer_boundary_nm"
    ]

    left_outer = axial_boundaries[
        "left_outer_boundary_nm"
    ]

    left_inner = axial_boundaries[
        "left_inner_boundary_nm"
    ]

    right_inner = axial_boundaries[
        "right_inner_boundary_nm"
    ]

    right_outer = axial_boundaries[
        "right_outer_boundary_nm"
    ]

    mouth_mask = (
        (
            (
                z_grid >= left_outer
            )
            & (
                z_grid < left_inner
            )
        )
        |
        (
            (
                z_grid > right_inner
            )
            & (
                z_grid <= right_outer
            )
        )
    ) & (
        r_grid < outer_r
    )

    lumen_mask = (
        z_grid >= left_inner
    ) & (
        z_grid <= right_inner
    ) & (
        r_grid < inner_r
    )

    interfacial_mask = (
        z_grid >= left_inner
    ) & (
        z_grid <= right_inner
    ) & (
        r_grid >= inner_r
    ) & (
        r_grid < outer_r
    )

    exterior_mask = ~(
        mouth_mask
        | lumen_mask
        | interfacial_mask
    )

    masks = {
        "lumen_core": lumen_mask,
        "interfacial_shell": (
            interfacial_mask
        ),
        "mouth_transitions": (
            mouth_mask
        ),
        "exterior": exterior_mask,
    }

    coverage = np.zeros_like(
        z_grid,
        dtype=np.int8,
    )

    for mask in masks.values():
        coverage += mask.astype(
            np.int8
        )

    if not np.all(
        coverage == 1
    ):
        raise RuntimeError(
            "Profile-guided masks do not "
            "form an exact partition"
        )

    return masks


def summarize_regions(
    masks: dict[str, np.ndarray],
    cell_volumes: np.ndarray,
    integrated_counts: np.ndarray,
    total_volume: float,
    total_count: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    reconstructed_volume = 0.0
    reconstructed_count = 0.0

    for region_name in REGION_NAMES:
        mask = masks[
            region_name
        ]

        volume = float(
            np.sum(
                cell_volumes[
                    mask
                ]
            )
        )

        count = float(
            np.sum(
                integrated_counts[
                    mask
                ]
            )
        )

        reconstructed_volume += volume
        reconstructed_count += count

        rows.append(
            {
                "region": region_name,
                "display_name": (
                    REGION_LABELS[
                        region_name
                    ]
                ),
                "cell_count": int(
                    np.sum(mask)
                ),
                "volume_nm^3": volume,
                "volume_fraction": (
                    volume
                    / total_volume
                ),
                "average_water_count": (
                    count
                ),
                "water_count_fraction": (
                    count
                    / total_count
                ),
                "volume_weighted_density_nm^-3": (
                    count / volume
                ),
            }
        )

    count_error = abs(
        reconstructed_count
        - total_count
    )

    volume_error = abs(
        reconstructed_volume
        - total_volume
    )

    if count_error > (
        CONSERVATION_TOLERANCE
    ):
        raise RuntimeError(
            "Water-count conservation failed: "
            f"{count_error:.12e}"
        )

    if volume_error > (
        CONSERVATION_TOLERANCE
    ):
        raise RuntimeError(
            "Volume conservation failed: "
            f"{volume_error:.12e}"
        )

    return rows


def region_code_matrix(
    masks: dict[str, np.ndarray],
) -> np.ndarray:
    shape = next(
        iter(
            masks.values()
        )
    ).shape

    codes = np.zeros(
        shape,
        dtype=np.int8,
    )

    for region_name in REGION_NAMES:
        codes[
            masks[
                region_name
            ]
        ] = REGION_CODES[
            region_name
        ]

    if np.any(
        codes == 0
    ):
        raise RuntimeError(
            "Unclassified density cells remain"
        )

    return codes


def load_operational_summary(
) -> dict[str, dict[str, float]]:
    rows = read_csv(
        OPERATIONAL_SUMMARY_CSV
    )

    mapping = {
        "interior_core": "lumen_core",
        "interfacial_shell": (
            "interfacial_shell"
        ),
        "mouth_zones": (
            "mouth_transitions"
        ),
        "exterior": "exterior",
    }

    result: dict[
        str,
        dict[str, float]
    ] = {}

    for row in rows:
        original_name = row["region"]

        if original_name not in mapping:
            continue

        result[
            mapping[
                original_name
            ]
        ] = {
            "count": float(
                row[
                    "average_water_count"
                ]
            ),
            "fraction": float(
                row[
                    "water_count_fraction"
                ]
            ),
        }

    missing = [
        region
        for region in REGION_NAMES
        if region not in result
    ]

    if missing:
        raise RuntimeError(
            "Operational comparison is missing: "
            + ", ".join(missing)
        )

    return result


def build_figure(
    density_data: dict[str, np.ndarray | float],
    pyrenes: list[dict[str, float | str]],
    radial_profile_raw: np.ndarray,
    radial_profile_smooth: np.ndarray,
    axial_profile_raw: np.ndarray,
    axial_profile_smooth: np.ndarray,
    radial_boundaries: dict[str, float],
    axial_boundaries: dict[str, float],
    main_summary: list[dict[str, object]],
    operational_summary: dict[
        str,
        dict[str, float],
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

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 8.8,
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
        figsize=(13.5, 9.8),
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

    inner_r = radial_boundaries[
        "inner_boundary_nm"
    ]

    outer_r = radial_boundaries[
        "outer_boundary_nm"
    ]

    minimum_r = radial_boundaries[
        "minimum_radius_nm"
    ]

    for radius, style in (
        (inner_r, "--"),
        (minimum_r, "-"),
        (outer_r, "--"),
    ):
        axis_a.axhline(
            radius,
            color="white",
            linestyle=style,
            linewidth=1.2,
        )

    axial_lines = (
        axial_boundaries[
            "left_outer_boundary_nm"
        ],
        axial_boundaries[
            "left_inner_boundary_nm"
        ],
        axial_boundaries[
            "right_inner_boundary_nm"
        ],
        axial_boundaries[
            "right_outer_boundary_nm"
        ],
    )

    for z_value in axial_lines:
        axis_a.axvline(
            z_value,
            color="white",
            linestyle=":",
            linewidth=1.2,
        )

    for pyrene in pyrenes:
        axis_a.scatter(
            float(
                pyrene["z_nm"]
            ),
            float(
                pyrene["r_nm"]
            ),
            marker="D",
            s=28,
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
        "Profile-guided spatial boundaries"
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
        radial_profile_raw,
        linewidth=1.0,
        alpha=0.45,
        label="Raw profile",
    )

    axis_b.plot(
        radial_centers,
        radial_profile_smooth,
        linewidth=2.0,
        label="Smoothed profile",
    )

    axis_b.axvline(
        minimum_r,
        color="black",
        linestyle="-",
        label="Depletion minimum",
    )

    axis_b.axvline(
        inner_r,
        color="black",
        linestyle="--",
    )

    axis_b.axvline(
        outer_r,
        color="black",
        linestyle="--",
    )

    axis_b.axhline(
        radial_boundaries[
            "inner_target_density_nm^-3"
        ],
        linestyle=":",
        linewidth=1.0,
    )

    axis_b.axhline(
        radial_boundaries[
            "outer_target_density_nm^-3"
        ],
        linestyle=":",
        linewidth=1.0,
    )

    axis_b.axvspan(
        inner_r,
        outer_r,
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
        "Radial depletion-shell detection"
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
        axial_profile_raw,
        linewidth=1.0,
        alpha=0.45,
        label="Raw profile",
    )

    axis_c.plot(
        axial_centers,
        axial_profile_smooth,
        linewidth=2.0,
        label="Smoothed profile",
    )

    left_outer = axial_boundaries[
        "left_outer_boundary_nm"
    ]

    left_inner = axial_boundaries[
        "left_inner_boundary_nm"
    ]

    right_inner = axial_boundaries[
        "right_inner_boundary_nm"
    ]

    right_outer = axial_boundaries[
        "right_outer_boundary_nm"
    ]

    axis_c.axvspan(
        left_outer,
        left_inner,
        alpha=0.12,
        label="Mouth transitions",
    )

    axis_c.axvspan(
        right_inner,
        right_outer,
        alpha=0.12,
    )

    for boundary in (
        left_outer,
        left_inner,
        right_inner,
        right_outer,
    ):
        axis_c.axvline(
            boundary,
            color="black",
            linestyle=":",
            linewidth=1.0,
        )

    axis_c.set_xlabel(
        "Axial coordinate z (nm)"
    )

    axis_c.set_ylabel(
        r"Lumen-averaged density "
        r"(nm$^{-3}$)"
    )

    axis_c.set_title(
        "Axial mouth-transition detection"
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

    x_positions = np.arange(
        len(REGION_NAMES)
    )

    profile_counts = np.asarray(
        [
            float(
                next(
                    row[
                        "average_water_count"
                    ]
                    for row
                    in main_summary
                    if row["region"]
                    == region
                )
            )
            for region in REGION_NAMES
        ],
        dtype=np.float64,
    )

    operational_counts = np.asarray(
        [
            operational_summary[
                region
            ]["count"]
            for region in REGION_NAMES
        ],
        dtype=np.float64,
    )

    bar_width = 0.38

    axis_d.bar(
        x_positions
        - 0.5 * bar_width,
        operational_counts,
        width=bar_width,
        label="Operational widths",
    )

    axis_d.bar(
        x_positions
        + 0.5 * bar_width,
        profile_counts,
        width=bar_width,
        label="Profile-guided",
    )

    axis_d.set_xticks(
        x_positions
    )

    axis_d.set_xticklabels(
        [
            REGION_LABELS[
                region
            ]
            for region in REGION_NAMES
        ],
        rotation=18,
        ha="right",
    )

    axis_d.set_ylabel(
        "Average water molecules"
    )

    axis_d.set_title(
        "Operational versus profile-guided partition"
    )

    axis_d.grid(
        True,
        axis="y",
        alpha=0.22,
    )

    axis_d.legend(
        frameon=False,
    )

    maximum_count = float(
        max(
            np.max(profile_counts),
            np.max(operational_counts),
        )
    )

    axis_d.set_ylim(
        0.0,
        maximum_count * 1.16,
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
            "Profile-guided confined-water "
            "classification around the continuous "
            "HBN scaffold"
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

    density_data = load_density()
    hbn_geometry = load_hbn_geometry()
    pyrenes = load_pyrenes()
    operational_summary = (
        load_operational_summary()
    )

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

    segment_lower = hbn_geometry[
        "lower_boundary_nm"
    ]

    segment_upper = hbn_geometry[
        "upper_boundary_nm"
    ]

    wall_radius = hbn_geometry[
        "mean_radius_nm"
    ]

    central_profile_mask = (
        axial_centers
        >= (
            segment_lower
            + END_EXCLUSION_NM
        )
    ) & (
        axial_centers
        <= (
            segment_upper
            - END_EXCLUSION_NM
        )
    )

    radial_profile_raw = (
        axial_average_radial_profile(
            density,
            axial_edges,
            central_profile_mask,
        )
    )

    radial_spacing = float(
        np.median(
            np.diff(
                radial_centers
            )
        )
    )

    radial_profile_smooth = (
        gaussian_smooth(
            radial_profile_raw,
            radial_spacing,
            RADIAL_SMOOTHING_SIGMA_NM,
        )
    )

    radial_boundaries = (
        derive_radial_boundaries(
            radial_centers,
            radial_profile_smooth,
            wall_radius,
            MAIN_RADIAL_RECOVERY_FRACTION,
        )
    )

    lumen_radial_mask = (
        radial_centers
        < radial_boundaries[
            "inner_boundary_nm"
        ]
    )

    axial_profile_raw = (
        radial_average_axial_profile(
            density,
            radial_edges,
            lumen_radial_mask,
        )
    )

    axial_spacing = float(
        np.median(
            np.diff(
                axial_centers
            )
        )
    )

    axial_profile_smooth = (
        gaussian_smooth(
            axial_profile_raw,
            axial_spacing,
            AXIAL_SMOOTHING_SIGMA_NM,
        )
    )

    axial_boundaries = (
        derive_axial_boundaries(
            axial_centers,
            axial_profile_smooth,
            segment_lower,
            segment_upper,
            MAIN_AXIAL_INNER_FRACTION,
            MAIN_AXIAL_OUTER_FRACTION,
        )
    )

    main_masks = classify_regions(
        axial_centers,
        radial_centers,
        radial_boundaries,
        axial_boundaries,
    )

    main_summary = summarize_regions(
        main_masks,
        cell_volumes,
        integrated_counts,
        total_volume,
        total_count,
    )

    for row in main_summary:
        row[
            "radial_recovery_fraction"
        ] = (
            MAIN_RADIAL_RECOVERY_FRACTION
        )

        row[
            "axial_inner_fraction"
        ] = (
            MAIN_AXIAL_INNER_FRACTION
        )

        row[
            "axial_outer_fraction"
        ] = (
            MAIN_AXIAL_OUTER_FRACTION
        )

    write_csv(
        REGION_SUMMARY_CSV,
        main_summary,
    )

    boundary_rows = [
        {
            "boundary_type": "radial",
            **radial_boundaries,
            "HBN_mean_radius_nm": (
                wall_radius
            ),
            "HBN_p05_radius_nm": (
                hbn_geometry[
                    "p05_radius_nm"
                ]
            ),
            "HBN_p95_radius_nm": (
                hbn_geometry[
                    "p95_radius_nm"
                ]
            ),
        },
        {
            "boundary_type": "axial",
            **axial_boundaries,
            "HBN_lower_boundary_nm": (
                segment_lower
            ),
            "HBN_upper_boundary_nm": (
                segment_upper
            ),
        },
    ]

    write_csv(
        BOUNDARIES_CSV,
        boundary_rows,
    )

    radial_profile_rows = [
        {
            "radial_coordinate_nm": (
                radial_centers[index]
            ),
            "raw_density_nm^-3": (
                radial_profile_raw[index]
            ),
            "smoothed_density_nm^-3": (
                radial_profile_smooth[
                    index
                ]
            ),
        }
        for index
        in range(
            radial_centers.size
        )
    ]

    write_csv(
        RADIAL_PROFILE_CSV,
        radial_profile_rows,
    )

    axial_profile_rows = [
        {
            "axial_coordinate_nm": (
                axial_centers[index]
            ),
            "raw_density_nm^-3": (
                axial_profile_raw[index]
            ),
            "smoothed_density_nm^-3": (
                axial_profile_smooth[
                    index
                ]
            ),
        }
        for index
        in range(
            axial_centers.size
        )
    ]

    write_csv(
        AXIAL_PROFILE_CSV,
        axial_profile_rows,
    )

    sensitivity_rows: list[
        dict[str, object]
    ] = []

    for recovery_fraction in (
        RADIAL_RECOVERY_FRACTIONS
    ):
        current_radial = (
            derive_radial_boundaries(
                radial_centers,
                radial_profile_smooth,
                wall_radius,
                recovery_fraction,
            )
        )

        current_lumen_mask = (
            radial_centers
            < current_radial[
                "inner_boundary_nm"
            ]
        )

        current_axial_raw = (
            radial_average_axial_profile(
                density,
                radial_edges,
                current_lumen_mask,
            )
        )

        current_axial_smooth = (
            gaussian_smooth(
                current_axial_raw,
                axial_spacing,
                AXIAL_SMOOTHING_SIGMA_NM,
            )
        )

        for inner_fraction in (
            AXIAL_INNER_FRACTIONS
        ):
            outer_fraction = (
                1.0
                - inner_fraction
            )

            current_axial = (
                derive_axial_boundaries(
                    axial_centers,
                    current_axial_smooth,
                    segment_lower,
                    segment_upper,
                    inner_fraction,
                    outer_fraction,
                )
            )

            current_masks = classify_regions(
                axial_centers,
                radial_centers,
                current_radial,
                current_axial,
            )

            current_summary = (
                summarize_regions(
                    current_masks,
                    cell_volumes,
                    integrated_counts,
                    total_volume,
                    total_count,
                )
            )

            for row in current_summary:
                sensitivity_rows.append(
                    {
                        "radial_recovery_fraction": (
                            recovery_fraction
                        ),
                        "axial_inner_fraction": (
                            inner_fraction
                        ),
                        "axial_outer_fraction": (
                            outer_fraction
                        ),
                        "inner_radial_boundary_nm": (
                            current_radial[
                                "inner_boundary_nm"
                            ]
                        ),
                        "outer_radial_boundary_nm": (
                            current_radial[
                                "outer_boundary_nm"
                            ]
                        ),
                        "left_outer_boundary_nm": (
                            current_axial[
                                "left_outer_boundary_nm"
                            ]
                        ),
                        "left_inner_boundary_nm": (
                            current_axial[
                                "left_inner_boundary_nm"
                            ]
                        ),
                        "right_inner_boundary_nm": (
                            current_axial[
                                "right_inner_boundary_nm"
                            ]
                        ),
                        "right_outer_boundary_nm": (
                            current_axial[
                                "right_outer_boundary_nm"
                            ]
                        ),
                        **row,
                    }
                )

    write_csv(
        SENSITIVITY_CSV,
        sensitivity_rows,
    )

    main_lookup = {
        str(row["region"]): row
        for row in main_summary
    }

    comparison_rows = []

    for region_name in REGION_NAMES:
        operational = (
            operational_summary[
                region_name
            ]
        )

        profile = main_lookup[
            region_name
        ]

        profile_count = float(
            profile[
                "average_water_count"
            ]
        )

        profile_fraction = float(
            profile[
                "water_count_fraction"
            ]
        )

        comparison_rows.append(
            {
                "region": region_name,
                "display_name": (
                    REGION_LABELS[
                        region_name
                    ]
                ),
                "operational_count": (
                    operational["count"]
                ),
                "profile_guided_count": (
                    profile_count
                ),
                "count_difference": (
                    profile_count
                    - operational["count"]
                ),
                "relative_count_change": (
                    (
                        profile_count
                        - operational["count"]
                    )
                    / operational["count"]
                ),
                "operational_fraction": (
                    operational["fraction"]
                ),
                "profile_guided_fraction": (
                    profile_fraction
                ),
                "fraction_difference": (
                    profile_fraction
                    - operational["fraction"]
                ),
            }
        )

    write_csv(
        COMPARISON_CSV,
        comparison_rows,
    )

    region_codes = (
        region_code_matrix(
            main_masks
        )
    )

    np.savez_compressed(
        REGIONAL_NPZ,
        axial_edges_nm=axial_edges,
        radial_edges_nm=radial_edges,
        axial_centers_nm=(
            axial_centers
        ),
        radial_centers_nm=(
            radial_centers
        ),
        density_nm3=density,
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
                REGION_CODES[
                    region
                ]
                for region
                in REGION_NAMES
            ],
            dtype=np.int8,
        ),
        radial_profile_raw_nm3=(
            radial_profile_raw
        ),
        radial_profile_smooth_nm3=(
            radial_profile_smooth
        ),
        axial_profile_raw_nm3=(
            axial_profile_raw
        ),
        axial_profile_smooth_nm3=(
            axial_profile_smooth
        ),
        inner_radial_boundary_nm=np.asarray(
            radial_boundaries[
                "inner_boundary_nm"
            ]
        ),
        radial_minimum_nm=np.asarray(
            radial_boundaries[
                "minimum_radius_nm"
            ]
        ),
        outer_radial_boundary_nm=np.asarray(
            radial_boundaries[
                "outer_boundary_nm"
            ]
        ),
        left_outer_boundary_nm=np.asarray(
            axial_boundaries[
                "left_outer_boundary_nm"
            ]
        ),
        left_inner_boundary_nm=np.asarray(
            axial_boundaries[
                "left_inner_boundary_nm"
            ]
        ),
        right_inner_boundary_nm=np.asarray(
            axial_boundaries[
                "right_inner_boundary_nm"
            ]
        ),
        right_outer_boundary_nm=np.asarray(
            axial_boundaries[
                "right_outer_boundary_nm"
            ]
        ),
        total_integrated_water_count=np.asarray(
            total_count
        ),
    )

    build_figure(
        density_data,
        pyrenes,
        radial_profile_raw,
        radial_profile_smooth,
        axial_profile_raw,
        axial_profile_smooth,
        radial_boundaries,
        axial_boundaries,
        main_summary,
        operational_summary,
    )

    reconstructed_count = sum(
        float(
            row[
                "average_water_count"
            ]
        )
        for row in main_summary
    )

    conservation_error = abs(
        reconstructed_count
        - total_count
    )

    sensitivity_by_region: dict[
        str,
        list[float],
    ] = {
        region: []
        for region in REGION_NAMES
    }

    for row in sensitivity_rows:
        sensitivity_by_region[
            str(row["region"])
        ].append(
            float(
                row[
                    "average_water_count"
                ]
            )
        )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Profile-Guided "
            "Confined-Water Classification\n\n"
        )

        handle.write(
            "## Geometry and profile basis\n\n"
        )

        handle.write(
            "The HBN scaffold is one continuous "
            "segment. Regional boundaries were "
            "derived from the smoothed radial and "
            "axial water-density profiles rather "
            "than from fixed geometric widths.\n\n"
        )

        handle.write(
            f"- HBN axial extent: "
            f"{segment_lower:.6f} to "
            f"{segment_upper:.6f} nm.\n"
        )

        handle.write(
            f"- Mean HBN radius: "
            f"{wall_radius:.6f} nm.\n"
        )

        handle.write(
            f"- Radial depletion minimum: "
            f"{radial_boundaries['minimum_radius_nm']:.6f} nm.\n"
        )

        handle.write(
            f"- Inner radial boundary: "
            f"{radial_boundaries['inner_boundary_nm']:.6f} nm.\n"
        )

        handle.write(
            f"- Outer radial boundary: "
            f"{radial_boundaries['outer_boundary_nm']:.6f} nm.\n"
        )

        handle.write(
            f"- Left mouth transition: "
            f"{axial_boundaries['left_outer_boundary_nm']:.6f} to "
            f"{axial_boundaries['left_inner_boundary_nm']:.6f} nm.\n"
        )

        handle.write(
            f"- Right mouth transition: "
            f"{axial_boundaries['right_inner_boundary_nm']:.6f} to "
            f"{axial_boundaries['right_outer_boundary_nm']:.6f} nm.\n\n"
        )

        handle.write(
            "## Profile-guided populations\n\n"
        )

        for row in main_summary:
            handle.write(
                f"- {row['display_name']}: "
                f"{float(row['average_water_count']):.6f} "
                "waters "
                f"({100.0 * float(row['water_count_fraction']):.3f}%).\n"
            )

        handle.write(
            "\n## Conservation\n\n"
        )

        handle.write(
            f"- Canonical cylinder count: "
            f"{total_count:.9f}.\n"
        )

        handle.write(
            f"- Reconstructed count: "
            f"{reconstructed_count:.9f}.\n"
        )

        handle.write(
            f"- Absolute error: "
            f"{conservation_error:.3e}.\n\n"
        )

        handle.write(
            "## Boundary sensitivity\n\n"
        )

        handle.write(
            "The analysis was repeated for radial "
            "recovery fractions of 0.40, 0.50, and "
            "0.60 and axial inner fractions of "
            "0.75, 0.80, and 0.85.\n\n"
        )

        for region_name in REGION_NAMES:
            values = sensitivity_by_region[
                region_name
            ]

            central = float(
                main_lookup[
                    region_name
                ][
                    "average_water_count"
                ]
            )

            span = (
                max(values)
                - min(values)
            )

            handle.write(
                f"- {REGION_LABELS[region_name]}: "
                f"{min(values):.6f} to "
                f"{max(values):.6f} waters; "
                f"span relative to central value "
                f"{100.0 * span / central:.3f}%.\n"
            )

        handle.write(
            "\n## Interpretation boundary\n\n"
        )

        handle.write(
            "The resulting populations are "
            "profile-guided effective occupancies. "
            "They characterize time-averaged solvent "
            "organization around a frozen solute. "
            "They are not residence-time populations "
            "and do not establish scaffold thermal "
            "stability.\n"
        )

    with MANIFEST_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Profile-Guided "
            "Analysis Manifest\n\n"
        )

        for path in (
            BOUNDARIES_CSV,
            REGION_SUMMARY_CSV,
            SENSITIVITY_CSV,
            COMPARISON_CSV,
            RADIAL_PROFILE_CSV,
            AXIAL_PROFILE_CSV,
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
        BOUNDARIES_CSV,
        REGION_SUMMARY_CSV,
        SENSITIVITY_CSV,
        COMPARISON_CSV,
        RADIAL_PROFILE_CSV,
        AXIAL_PROFILE_CSV,
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

    missing_outputs = [
        path
        for path in required_outputs
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing_outputs:
        raise RuntimeError(
            "Missing or empty outputs:\n"
            + "\n".join(
                str(path)
                for path in missing_outputs
            )
        )

    log(
        "Day020 profile-guided confined-water "
        "analysis completed."
    )

    log(
        "HBN architecture: "
        "one continuous segment"
    )

    log(
        "Radial depletion minimum: "
        f"{radial_boundaries['minimum_radius_nm']:.6f} nm"
    )

    log(
        "Radial boundaries: "
        f"{radial_boundaries['inner_boundary_nm']:.6f} to "
        f"{radial_boundaries['outer_boundary_nm']:.6f} nm"
    )

    log(
        "Left mouth transition: "
        f"{axial_boundaries['left_outer_boundary_nm']:.6f} to "
        f"{axial_boundaries['left_inner_boundary_nm']:.6f} nm"
    )

    log(
        "Right mouth transition: "
        f"{axial_boundaries['right_inner_boundary_nm']:.6f} to "
        f"{axial_boundaries['right_outer_boundary_nm']:.6f} nm"
    )

    log(
        f"Canonical cylinder count: "
        f"{total_count:.6f}"
    )

    for row in main_summary:
        log(
            f"{row['display_name']}: "
            f"{float(row['average_water_count']):.6f} "
            f"({100.0 * float(row['water_count_fraction']):.3f}%)"
        )

    log(
        f"Regional conservation error: "
        f"{conservation_error:.3e}"
    )

    log(
        f"Sensitivity conditions: "
        f"{len(RADIAL_RECOVERY_FRACTIONS) * len(AXIAL_INNER_FRACTIONS)}"
    )

    for region_name in REGION_NAMES:
        values = sensitivity_by_region[
            region_name
        ]

        central = float(
            main_lookup[
                region_name
            ][
                "average_water_count"
            ]
        )

        relative_span = (
            max(values)
            - min(values)
        ) / central

        log(
            f"{REGION_LABELS[region_name]} "
            f"sensitivity span: "
            f"{100.0 * relative_span:.2f}%"
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
