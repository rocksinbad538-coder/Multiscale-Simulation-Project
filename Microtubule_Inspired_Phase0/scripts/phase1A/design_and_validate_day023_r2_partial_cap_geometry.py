#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

R0_REFERENCE_ROOT = (
    DAY023_ROOT
    / "01_r0_t0_reference"
)

R1_PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

R1_50PS_ROOT = (
    DAY023_ROOT
    / "11_r1_frozen_solute_nvt_20_to_50ps"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "12_r2_partial_cap_geometry_design"
)

SELECTED_ROOT = (
    OUTPUT_ROOT
    / "selected"
)

R0_T0_GRO = (
    R0_REFERENCE_ROOT
    / "r0_accepted_t0_hydrated_system.gro"
)

R1_CAPS_GRO = (
    R1_PROTOTYPE_ROOT
    / "r1_selected_steric_caps_only.gro"
)

R1_PROTOTYPE_JSON = (
    R1_PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

R1_50PS_SUMMARY = (
    R1_50PS_ROOT
    / "r1_frozen_solute_nvt_50ps_summary.csv"
)

CANDIDATE_SCAN_CSV = (
    OUTPUT_ROOT
    / "r2_partial_cap_candidate_scan.csv"
)

SELECTED_CAP_COORDINATES_CSV = (
    SELECTED_ROOT
    / "r2_selected_partial_cap_coordinates.csv"
)

SELECTED_REMOVED_WATERS_CSV = (
    SELECTED_ROOT
    / "r2_selected_removed_water_molecules.csv"
)

SELECTED_CAPS_GRO = (
    SELECTED_ROOT
    / "r2_selected_partial_caps_only.gro"
)

SELECTED_SYSTEM_GRO = (
    SELECTED_ROOT
    / "r2_selected_partial_cap_geometry_only.gro"
)

SELECTED_DEFINITION_JSON = (
    SELECTED_ROOT
    / "r2_selected_partial_cap_definition.json"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_partial_cap_geometry_summary.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r2_partial_cap_geometry_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_PARTIAL_CAP_GEOMETRY_DESIGN_DAY023.md"
)

HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS

R0_WATERS = 16634
WATER_SITES = 4

R1_CAP_BEADS_PER_END = 163
R1_TOTAL_CAP_BEADS = 326

EXPECTED_R0_ATOMS = (
    SOLUTE_ATOMS
    + R0_WATERS * WATER_SITES
)

INITIAL_LUMEN_WATERS = 428

CAP_WATER_5KBT_DISTANCE_NM = 0.17
WATER_OVERLAP_CUTOFF_NM = 0.22

TARGET_EFFECTIVE_APERTURE_RADIUS_NM = 0.30

CANDIDATE_REMOVAL_RADII_NM = [
    round(
        0.20 + 0.05 * index,
        2,
    )
    for index in range(10)
]

MIN_CAP_BEADS_PER_END = 80

MIN_EFFECTIVE_APERTURE_RADIUS_NM = 0.20
MAX_EFFECTIVE_APERTURE_RADIUS_NM = 0.45

MIN_CONSERVATIVE_APERTURE_RADIUS_NM = 0.15
MAX_OPEN_AREA_FRACTION = 0.16

MIN_LUMEN_RETENTION_FRACTION = 0.98

MIN_CAP_HBN_DISTANCE_NM = 0.19
MIN_CAP_PYR_DISTANCE_NM = 0.90
MIN_CAP_OW_DISTANCE_NM = 0.2195


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require_file(path: Path) -> None:
    if (
        not path.exists()
        or not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; "
            f"found {len(rows)}"
        )

    return rows[0]


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

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
                    key: row.get(
                        key,
                        "",
                    )
                    for key in fields
                }
            )


def parse_box(
    values: list[float],
) -> np.ndarray:
    if len(values) != 3:
        raise RuntimeError(
            "This workflow requires an orthorhombic box."
        )

    box = np.asarray(
        values,
        dtype=float,
    )

    if (
        not np.all(
            np.isfinite(box)
        )
        or np.any(box <= 0.0)
    ):
        raise RuntimeError(
            "Invalid simulation box."
        )

    return box


def read_gro(
    path: Path,
) -> tuple[
    str,
    list[dict[str, Any]],
    np.ndarray,
]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    title = lines[0]
    atom_count = int(
        lines[1].strip()
    )

    if len(lines) < atom_count + 3:
        raise RuntimeError(
            f"Incomplete GRO file: {path}"
        )

    atoms: list[dict[str, Any]] = []

    for index, line in enumerate(
        lines[
            2:
            2 + atom_count
        ]
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed atom line "
                f"{index + 1} in {path}"
            )

        atom = {
            "index": index,
            "resid": int(
                line[0:5]
            ),
            "resname": line[
                5:10
            ].strip(),
            "atomname": line[
                10:15
            ].strip(),
            "atomnum": int(
                line[15:20]
            ),
            "position": np.asarray(
                [
                    float(
                        line[20:28]
                    ),
                    float(
                        line[28:36]
                    ),
                    float(
                        line[36:44]
                    ),
                ],
                dtype=float,
            ),
        }

        if len(line) >= 68:
            try:
                atom["velocity"] = np.asarray(
                    [
                        float(
                            line[44:52]
                        ),
                        float(
                            line[52:60]
                        ),
                        float(
                            line[60:68]
                        ),
                    ],
                    dtype=float,
                )
            except ValueError:
                atom["velocity"] = None
        else:
            atom["velocity"] = None

        atoms.append(atom)

    box = parse_box(
        [
            float(value)
            for value in lines[
                2 + atom_count
            ].split()
        ]
    )

    return (
        title,
        atoms,
        box,
    )


def write_gro(
    path: Path,
    title: str,
    atoms: list[dict[str, Any]],
    box: np.ndarray,
) -> None:
    lines = [
        title,
        f"{len(atoms):5d}",
    ]

    for atom_index, atom in enumerate(
        atoms,
        start=1,
    ):
        resid = int(
            atom["resid"]
        ) % 100000

        atomnum = atom_index % 100000

        if atomnum == 0:
            atomnum = 99999

        position = np.asarray(
            atom["position"],
            dtype=float,
        )

        line = (
            f"{resid:5d}"
            f"{str(atom['resname'])[:5]:<5}"
            f"{str(atom['atomname'])[:5]:>5}"
            f"{atomnum:5d}"
            f"{position[0]:8.3f}"
            f"{position[1]:8.3f}"
            f"{position[2]:8.3f}"
        )

        velocity = atom.get(
            "velocity"
        )

        if velocity is not None:
            velocity_array = np.asarray(
                velocity,
                dtype=float,
            )

            if np.all(
                np.isfinite(
                    velocity_array
                )
            ):
                line += (
                    f"{velocity_array[0]:8.4f}"
                    f"{velocity_array[1]:8.4f}"
                    f"{velocity_array[2]:8.4f}"
                )

        lines.append(line)

    lines.append(
        " ".join(
            f"{value:10.5f}"
            for value in box
        )
    )

    path.write_text(
        "\n".join(lines)
        + "\n",
        encoding="utf-8",
    )


def atom_positions(
    atoms: list[dict[str, Any]],
) -> np.ndarray:
    return np.asarray(
        [
            atom["position"]
            for atom in atoms
        ],
        dtype=float,
    )


def minimum_image(
    displacement: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    return (
        displacement
        - box
        * np.round(
            displacement
            / box
        )
    )


def pair_minimum_distance(
    first: np.ndarray,
    second: np.ndarray,
    box: np.ndarray,
    *,
    chunk_size: int = 512,
) -> float:
    if (
        len(first) == 0
        or len(second) == 0
    ):
        return math.inf

    minimum = math.inf

    for start in range(
        0,
        len(first),
        chunk_size,
    ):
        stop = min(
            start + chunk_size,
            len(first),
        )

        displacement = (
            first[
                start:stop,
                None,
                :,
            ]
            - second[
                None,
                :,
                :,
            ]
        )

        displacement = minimum_image(
            displacement,
            box,
        )

        distances = np.linalg.norm(
            displacement,
            axis=2,
        )

        minimum = min(
            minimum,
            float(
                np.min(distances)
            ),
        )

    return minimum


def nearest_distances(
    points: np.ndarray,
    reference: np.ndarray,
    box: np.ndarray,
    *,
    chunk_size: int = 512,
) -> np.ndarray:
    if len(reference) == 0:
        raise RuntimeError(
            "No reference atoms were supplied."
        )

    output = np.empty(
        len(points),
        dtype=float,
    )

    for start in range(
        0,
        len(points),
        chunk_size,
    ):
        stop = min(
            start + chunk_size,
            len(points),
        )

        displacement = (
            points[
                start:stop,
                None,
                :,
            ]
            - reference[
                None,
                :,
                :,
            ]
        )

        displacement = minimum_image(
            displacement,
            box,
        )

        distances = np.linalg.norm(
            displacement,
            axis=2,
        )

        output[
            start:stop
        ] = np.min(
            distances,
            axis=1,
        )

    return output


def water_oxygen_indices(
    atoms: list[dict[str, Any]],
) -> np.ndarray:
    indices = []

    for water_index in range(
        R0_WATERS
    ):
        start = (
            SOLUTE_ATOMS
            + water_index
            * WATER_SITES
        )

        chunk = atoms[
            start:
            start + WATER_SITES
        ]

        matches = [
            start + local_index
            for local_index, atom
            in enumerate(chunk)
            if atom["atomname"] == "OW"
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Could not identify exactly one OW "
                f"in water molecule {water_index}."
            )

        indices.append(
            matches[0]
        )

    return np.asarray(
        indices,
        dtype=int,
    )


def lumen_mask(
    water_positions: np.ndarray,
    box: np.ndarray,
    prototype: dict[str, Any],
) -> np.ndarray:
    center = np.asarray(
        prototype[
            "tube_center_wrapped_nm"
        ],
        dtype=float,
    )

    axis = np.asarray(
        prototype[
            "tube_axis"
        ],
        dtype=float,
    )

    axis /= np.linalg.norm(axis)

    axial_low = float(
        prototype[
            "axial_low_nm"
        ]
    )

    axial_high = float(
        prototype[
            "axial_high_nm"
        ]
    )

    accessible_radius = float(
        prototype[
            "accessible_radius_nm"
        ]
    )

    relative_coordinates = minimum_image(
        water_positions - center,
        box,
    )

    axial = (
        relative_coordinates
        @ axis
    )

    perpendicular = (
        relative_coordinates
        - np.outer(
            axial,
            axis,
        )
    )

    radial = np.linalg.norm(
        perpendicular,
        axis=1,
    )

    return (
        (axial >= axial_low)
        & (axial <= axial_high)
        & (radial <= accessible_radius)
    )


def cap_geometry(
    cap_positions: np.ndarray,
    box: np.ndarray,
    prototype: dict[str, Any],
) -> dict[str, np.ndarray]:
    center = np.asarray(
        prototype[
            "tube_center_wrapped_nm"
        ],
        dtype=float,
    )

    axis = np.asarray(
        prototype[
            "tube_axis"
        ],
        dtype=float,
    )

    axis /= np.linalg.norm(axis)

    relative_coordinates = minimum_image(
        cap_positions - center,
        box,
    )

    axial = (
        relative_coordinates
        @ axis
    )

    perpendicular = (
        relative_coordinates
        - np.outer(
            axial,
            axis,
        )
    )

    radial = np.linalg.norm(
        perpendicular,
        axis=1,
    )

    axial_midpoint = float(
        np.median(axial)
    )

    lower_mask = (
        axial < axial_midpoint
    )

    upper_mask = (
        axial > axial_midpoint
    )

    if (
        np.count_nonzero(lower_mask)
        != R1_CAP_BEADS_PER_END
        or np.count_nonzero(upper_mask)
        != R1_CAP_BEADS_PER_END
    ):
        raise RuntimeError(
            "Could not partition the R1 cap lattice into "
            "163 lower and 163 upper beads."
        )

    return {
        "axis": axis,
        "axial": axial,
        "radial": radial,
        "lower_mask": lower_mask,
        "upper_mask": upper_mask,
    }


def build_candidate(
    *,
    removal_radius_nm: float,
    cap_atoms: list[dict[str, Any]],
    cap_positions: np.ndarray,
    cap_data: dict[str, np.ndarray],
    water_oxygen_positions: np.ndarray,
    initial_lumen_mask: np.ndarray,
    hbn_positions: np.ndarray,
    pyr_positions: np.ndarray,
    box: np.ndarray,
    accessible_radius_nm: float,
) -> dict[str, Any]:
    radial = cap_data[
        "radial"
    ]

    lower_mask = cap_data[
        "lower_mask"
    ]

    upper_mask = cap_data[
        "upper_mask"
    ]

    retained_mask = (
        radial
        >= (
            removal_radius_nm
            - 1.0e-12
        )
    )

    retained_lower_mask = (
        retained_mask
        & lower_mask
    )

    retained_upper_mask = (
        retained_mask
        & upper_mask
    )

    retained_lower_indices = np.flatnonzero(
        retained_lower_mask
    )

    retained_upper_indices = np.flatnonzero(
        retained_upper_mask
    )

    retained_indices = np.concatenate(
        (
            retained_lower_indices,
            retained_upper_indices,
        )
    )

    retained_cap_positions = (
        cap_positions[
            retained_indices
        ]
    )

    lower_count = len(
        retained_lower_indices
    )

    upper_count = len(
        retained_upper_indices
    )

    if (
        lower_count == 0
        or upper_count == 0
    ):
        raise RuntimeError(
            "Candidate removed an entire cap."
        )

    lower_nearest_radial = float(
        np.min(
            radial[
                retained_lower_indices
            ]
        )
    )

    upper_nearest_radial = float(
        np.min(
            radial[
                retained_upper_indices
            ]
        )
    )

    limiting_nearest_radial = min(
        lower_nearest_radial,
        upper_nearest_radial,
    )

    effective_aperture_radius = max(
        0.0,
        limiting_nearest_radial
        - CAP_WATER_5KBT_DISTANCE_NM,
    )

    conservative_aperture_radius = max(
        0.0,
        limiting_nearest_radial
        - WATER_OVERLAP_CUTOFF_NM,
    )

    open_area_fraction = (
        effective_aperture_radius
        / accessible_radius_nm
    ) ** 2

    water_cap_distances = nearest_distances(
        water_oxygen_positions,
        retained_cap_positions,
        box,
    )

    remove_water_mask = (
        water_cap_distances
        < WATER_OVERLAP_CUTOFF_NM
    )

    retained_water_mask = (
        ~remove_water_mask
    )

    removed_water_count = int(
        np.count_nonzero(
            remove_water_mask
        )
    )

    retained_water_count = int(
        np.count_nonzero(
            retained_water_mask
        )
    )

    removed_lumen_water_count = int(
        np.count_nonzero(
            remove_water_mask
            & initial_lumen_mask
        )
    )

    retained_lumen_water_count = int(
        np.count_nonzero(
            retained_water_mask
            & initial_lumen_mask
        )
    )

    lumen_retention_fraction = (
        retained_lumen_water_count
        / INITIAL_LUMEN_WATERS
    )

    retained_water_positions = (
        water_oxygen_positions[
            retained_water_mask
        ]
    )

    minimum_cap_ow_distance = (
        pair_minimum_distance(
            retained_cap_positions,
            retained_water_positions,
            box,
        )
    )

    minimum_cap_hbn_distance = (
        pair_minimum_distance(
            retained_cap_positions,
            hbn_positions,
            box,
        )
    )

    minimum_cap_pyr_distance = (
        pair_minimum_distance(
            retained_cap_positions,
            pyr_positions,
            box,
        )
    )

    candidate_gates = {
        "symmetric_cap_counts": (
            lower_count
            == upper_count
        ),
        "minimum_cap_beads_per_end": (
            lower_count
            >= MIN_CAP_BEADS_PER_END
            and upper_count
            >= MIN_CAP_BEADS_PER_END
        ),
        "effective_aperture_radius_in_target_range": (
            MIN_EFFECTIVE_APERTURE_RADIUS_NM
            <= effective_aperture_radius
            <= MAX_EFFECTIVE_APERTURE_RADIUS_NM
        ),
        "conservative_aperture_is_open": (
            conservative_aperture_radius
            >= MIN_CONSERVATIVE_APERTURE_RADIUS_NM
        ),
        "open_area_fraction_is_limited": (
            open_area_fraction
            <= MAX_OPEN_AREA_FRACTION
        ),
        "lumen_water_retention_is_at_least_98_percent": (
            lumen_retention_fraction
            >= MIN_LUMEN_RETENTION_FRACTION
        ),
        "CAP_HBN_distance_is_valid": (
            minimum_cap_hbn_distance
            >= MIN_CAP_HBN_DISTANCE_NM
        ),
        "CAP_PYR_distance_is_valid": (
            minimum_cap_pyr_distance
            >= MIN_CAP_PYR_DISTANCE_NM
        ),
        "CAP_OW_distance_is_valid": (
            minimum_cap_ow_distance
            >= MIN_CAP_OW_DISTANCE_NM
        ),
    }

    failed_gates = [
        name
        for name, passed
        in candidate_gates.items()
        if not passed
    ]

    candidate_pass = (
        len(failed_gates) == 0
    )

    candidate_id = (
        "remove_r"
        + f"{removal_radius_nm:.2f}"
        .replace(
            ".",
            "p",
        )
    )

    score = (
        abs(
            effective_aperture_radius
            - TARGET_EFFECTIVE_APERTURE_RADIUS_NM
        ),
        -min(
            lower_count,
            upper_count,
        ),
        removed_water_count,
    )

    return {
        "candidate_id": candidate_id,
        "removal_radius_nm": (
            removal_radius_nm
        ),
        "retained_lower_indices": (
            retained_lower_indices
        ),
        "retained_upper_indices": (
            retained_upper_indices
        ),
        "retained_indices": (
            retained_indices
        ),
        "retained_water_mask": (
            retained_water_mask
        ),
        "remove_water_mask": (
            remove_water_mask
        ),
        "water_cap_distances": (
            water_cap_distances
        ),
        "cap_beads_lower": (
            lower_count
        ),
        "cap_beads_upper": (
            upper_count
        ),
        "cap_beads_total": (
            lower_count + upper_count
        ),
        "removed_cap_beads_total": (
            R1_TOTAL_CAP_BEADS
            - lower_count
            - upper_count
        ),
        "nearest_lower_cap_radial_nm": (
            lower_nearest_radial
        ),
        "nearest_upper_cap_radial_nm": (
            upper_nearest_radial
        ),
        "limiting_nearest_cap_radial_nm": (
            limiting_nearest_radial
        ),
        "effective_aperture_radius_5kBT_nm": (
            effective_aperture_radius
        ),
        "effective_aperture_diameter_5kBT_nm": (
            2.0
            * effective_aperture_radius
        ),
        "conservative_aperture_radius_nm": (
            conservative_aperture_radius
        ),
        "conservative_aperture_diameter_nm": (
            2.0
            * conservative_aperture_radius
        ),
        "open_area_fraction": (
            open_area_fraction
        ),
        "removed_water_molecules": (
            removed_water_count
        ),
        "retained_water_molecules": (
            retained_water_count
        ),
        "removed_lumen_water_molecules": (
            removed_lumen_water_count
        ),
        "retained_lumen_water_molecules": (
            retained_lumen_water_count
        ),
        "lumen_water_retention_fraction": (
            lumen_retention_fraction
        ),
        "minimum_CAP_HBN_distance_nm": (
            minimum_cap_hbn_distance
        ),
        "minimum_CAP_PYR_distance_nm": (
            minimum_cap_pyr_distance
        ),
        "minimum_CAP_OW_distance_nm": (
            minimum_cap_ow_distance
        ),
        "candidate_pass": (
            candidate_pass
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "selection_score": score,
    }


def candidate_csv_row(
    candidate: dict[str, Any],
) -> dict[str, Any]:
    excluded = {
        "retained_lower_indices",
        "retained_upper_indices",
        "retained_indices",
        "retained_water_mask",
        "remove_water_mask",
        "water_cap_distances",
        "selection_score",
    }

    return {
        key: value
        for key, value
        in candidate.items()
        if key not in excluded
    }


def clone_atom(
    atom: dict[str, Any],
) -> dict[str, Any]:
    return {
        "resid": int(
            atom["resid"]
        ),
        "resname": str(
            atom["resname"]
        ),
        "atomname": str(
            atom["atomname"]
        ),
        "position": np.asarray(
            atom["position"],
            dtype=float,
        ).copy(),
        "velocity": (
            None
            if atom.get(
                "velocity"
            ) is None
            else np.asarray(
                atom["velocity"],
                dtype=float,
            ).copy()
        ),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    SELECTED_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        R0_T0_GRO,
        R1_CAPS_GRO,
        R1_PROTOTYPE_JSON,
        R1_50PS_SUMMARY,
    ):
        require_file(required)

    r1_summary = read_single_csv_row(
        R1_50PS_SUMMARY
    )

    if (
        r1_summary.get(
            "decision"
        )
        !=
        "R1_FROZEN_SOLUTE_50PS_POSITIVE_CONTROL_VALIDATED"
    ):
        raise RuntimeError(
            "R1 has not been validated as the "
            "50 ps positive control."
        )

    if not parse_bool(
        r1_summary.get(
            "R2_static_design_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "R1 did not authorize the R2 static design."
        )

    prototype = json.loads(
        R1_PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    (
        r0_title,
        r0_atoms,
        r0_box,
    ) = read_gro(
        R0_T0_GRO
    )

    (
        caps_title,
        cap_atoms,
        cap_box,
    ) = read_gro(
        R1_CAPS_GRO
    )

    if len(r0_atoms) != EXPECTED_R0_ATOMS:
        raise RuntimeError(
            "Unexpected R0 atom count: "
            f"{len(r0_atoms)}/"
            f"{EXPECTED_R0_ATOMS}"
        )

    if len(cap_atoms) != R1_TOTAL_CAP_BEADS:
        raise RuntimeError(
            "Unexpected R1 cap-bead count: "
            f"{len(cap_atoms)}/"
            f"{R1_TOTAL_CAP_BEADS}"
        )

    if not np.allclose(
        r0_box,
        cap_box,
        rtol=0.0,
        atol=1.0e-6,
    ):
        raise RuntimeError(
            "R0 and R1 cap boxes do not match."
        )

    r0_positions = atom_positions(
        r0_atoms
    )

    cap_positions = atom_positions(
        cap_atoms
    )

    hbn_positions = (
        r0_positions[
            :HBN_ATOMS
        ]
    )

    pyr_positions = (
        r0_positions[
            HBN_ATOMS:
            SOLUTE_ATOMS
        ]
    )

    oxygen_indices = water_oxygen_indices(
        r0_atoms
    )

    water_oxygen_positions = (
        r0_positions[
            oxygen_indices
        ]
    )

    initial_lumen_mask = lumen_mask(
        water_oxygen_positions,
        r0_box,
        prototype,
    )

    initial_lumen_count = int(
        np.count_nonzero(
            initial_lumen_mask
        )
    )

    if (
        initial_lumen_count
        != INITIAL_LUMEN_WATERS
    ):
        raise RuntimeError(
            "The authoritative R0 lumen occupancy "
            "is not 428 waters: "
            f"{initial_lumen_count}"
        )

    cap_data = cap_geometry(
        cap_positions,
        r0_box,
        prototype,
    )

    accessible_radius_nm = float(
        prototype[
            "accessible_radius_nm"
        ]
    )

    candidates = []

    for removal_radius_nm in (
        CANDIDATE_REMOVAL_RADII_NM
    ):
        candidate = build_candidate(
            removal_radius_nm=(
                removal_radius_nm
            ),
            cap_atoms=cap_atoms,
            cap_positions=cap_positions,
            cap_data=cap_data,
            water_oxygen_positions=(
                water_oxygen_positions
            ),
            initial_lumen_mask=(
                initial_lumen_mask
            ),
            hbn_positions=(
                hbn_positions
            ),
            pyr_positions=(
                pyr_positions
            ),
            box=r0_box,
            accessible_radius_nm=(
                accessible_radius_nm
            ),
        )

        candidates.append(candidate)

    write_csv(
        CANDIDATE_SCAN_CSV,
        [
            candidate_csv_row(
                candidate
            )
            for candidate in candidates
        ],
    )

    passing_candidates = [
        candidate
        for candidate in candidates
        if candidate[
            "candidate_pass"
        ]
    ]

    if not passing_candidates:
        raise RuntimeError(
            "No R2 partial-cap geometry candidate "
            "passed the static gates."
        )

    selected = min(
        passing_candidates,
        key=lambda candidate: (
            candidate[
                "selection_score"
            ]
        ),
    )

    retained_water_mask = selected[
        "retained_water_mask"
    ]

    remove_water_mask = selected[
        "remove_water_mask"
    ]

    retained_lower_indices = selected[
        "retained_lower_indices"
    ]

    retained_upper_indices = selected[
        "retained_upper_indices"
    ]

    selected_atoms: list[
        dict[str, Any]
    ] = []

    for atom in r0_atoms[
        :SOLUTE_ATOMS
    ]:
        selected_atoms.append(
            clone_atom(atom)
        )

    removed_water_rows = []

    retained_water_counter = 0

    for water_index in range(
        R0_WATERS
    ):
        start = (
            SOLUTE_ATOMS
            + water_index
            * WATER_SITES
        )

        stop = (
            start
            + WATER_SITES
        )

        if remove_water_mask[
            water_index
        ]:
            oxygen_position = (
                water_oxygen_positions[
                    water_index
                ]
            )

            removed_water_rows.append(
                {
                    "source_water_index_zero_based": (
                        water_index
                    ),
                    "source_water_index_one_based": (
                        water_index + 1
                    ),
                    "source_oxygen_atom_one_based": (
                        int(
                            oxygen_indices[
                                water_index
                            ]
                        )
                        + 1
                    ),
                    "oxygen_x_nm": (
                        oxygen_position[0]
                    ),
                    "oxygen_y_nm": (
                        oxygen_position[1]
                    ),
                    "oxygen_z_nm": (
                        oxygen_position[2]
                    ),
                    "nearest_CAP_distance_nm": (
                        selected[
                            "water_cap_distances"
                        ][
                            water_index
                        ]
                    ),
                    "initially_luminal": bool(
                        initial_lumen_mask[
                            water_index
                        ]
                    ),
                }
            )

            continue

        retained_water_counter += 1

        for atom in r0_atoms[
            start:stop
        ]:
            selected_atoms.append(
                clone_atom(atom)
            )

    if (
        retained_water_counter
        != selected[
            "retained_water_molecules"
        ]
    ):
        raise RuntimeError(
            "Retained-water accounting failed."
        )

    last_water_resid = int(
        selected_atoms[
            -1
        ][
            "resid"
        ]
    )

    lower_resid = (
        last_water_resid
        + 1
    )

    upper_resid = (
        last_water_resid
        + 2
    )

    selected_cap_coordinate_rows = []

    selected_caps_atoms = []

    for end_name, indices, resid, resname in (
        (
            "lower",
            retained_lower_indices,
            lower_resid,
            "CAPL",
        ),
        (
            "upper",
            retained_upper_indices,
            upper_resid,
            "CAPU",
        ),
    ):
        for local_index, source_index in enumerate(
            indices,
            start=1,
        ):
            source_atom = cap_atoms[
                int(
                    source_index
                )
            ]

            position = np.asarray(
                source_atom[
                    "position"
                ],
                dtype=float,
            )

            cap_atom = {
                "resid": resid,
                "resname": resname,
                "atomname": "CAP",
                "position": position.copy(),
                "velocity": None,
            }

            selected_caps_atoms.append(
                cap_atom
            )

            selected_atoms.append(
                clone_atom(
                    cap_atom
                )
            )

            selected_cap_coordinate_rows.append(
                {
                    "end": end_name,
                    "end_local_index": (
                        local_index
                    ),
                    "source_R1_cap_index_zero_based": (
                        int(
                            source_index
                        )
                    ),
                    "source_R1_cap_index_one_based": (
                        int(
                            source_index
                        )
                        + 1
                    ),
                    "x_nm": position[0],
                    "y_nm": position[1],
                    "z_nm": position[2],
                    "radial_distance_from_axis_nm": (
                        cap_data[
                            "radial"
                        ][
                            int(
                                source_index
                            )
                        ]
                    ),
                    "axial_coordinate_nm": (
                        cap_data[
                            "axial"
                        ][
                            int(
                                source_index
                            )
                        ]
                    ),
                }
            )

    expected_selected_atoms = (
        SOLUTE_ATOMS
        + selected[
            "retained_water_molecules"
        ]
        * WATER_SITES
        + selected[
            "cap_beads_total"
        ]
    )

    if len(selected_atoms) != expected_selected_atoms:
        raise RuntimeError(
            "Selected R2 atom accounting failed: "
            f"{len(selected_atoms)}/"
            f"{expected_selected_atoms}"
        )

    write_gro(
        SELECTED_CAPS_GRO,
        (
            "R2 selected symmetric partial caps; "
            "geometry only"
        ),
        selected_caps_atoms,
        r0_box,
    )

    write_gro(
        SELECTED_SYSTEM_GRO,
        (
            "R2 hydrated partial-cap geometry; "
            "no topology or MD authorization"
        ),
        selected_atoms,
        r0_box,
    )

    write_csv(
        SELECTED_CAP_COORDINATES_CSV,
        selected_cap_coordinate_rows,
    )

    if removed_water_rows:
        write_csv(
            SELECTED_REMOVED_WATERS_CSV,
            removed_water_rows,
        )
    else:
        SELECTED_REMOVED_WATERS_CSV.write_text(
            (
                "source_water_index_zero_based,"
                "source_water_index_one_based,"
                "source_oxygen_atom_one_based,"
                "oxygen_x_nm,oxygen_y_nm,oxygen_z_nm,"
                "nearest_CAP_distance_nm,"
                "initially_luminal\n"
            ),
            encoding="utf-8",
        )

    selected_definition = {
        "decision": (
            "R2_PARTIAL_CAP_GEOMETRY_SELECTED"
        ),
        "architecture_label": (
            "R2 symmetric partial caps with "
            "central axial apertures"
        ),
        "model_class": (
            "neutral frozen steric screening design"
        ),
        "chemically_realizable_final_model": False,
        "source_R0_geometry": relative(
            R0_T0_GRO
        ),
        "source_R1_cap_lattice": relative(
            R1_CAPS_GRO
        ),
        "source_R1_validation_summary": relative(
            R1_50PS_SUMMARY
        ),
        "candidate_id": selected[
            "candidate_id"
        ],
        "selection_target_effective_aperture_radius_nm": (
            TARGET_EFFECTIVE_APERTURE_RADIUS_NM
        ),
        "removal_radius_nm": selected[
            "removal_radius_nm"
        ],
        "cap_water_5kBT_distance_nm": (
            CAP_WATER_5KBT_DISTANCE_NM
        ),
        "water_overlap_cutoff_nm": (
            WATER_OVERLAP_CUTOFF_NM
        ),
        "cap_beads_lower": selected[
            "cap_beads_lower"
        ],
        "cap_beads_upper": selected[
            "cap_beads_upper"
        ],
        "cap_beads_total": selected[
            "cap_beads_total"
        ],
        "removed_cap_beads_total": selected[
            "removed_cap_beads_total"
        ],
        "effective_aperture_radius_5kBT_nm": (
            selected[
                "effective_aperture_radius_5kBT_nm"
            ]
        ),
        "effective_aperture_diameter_5kBT_nm": (
            selected[
                "effective_aperture_diameter_5kBT_nm"
            ]
        ),
        "conservative_aperture_radius_nm": (
            selected[
                "conservative_aperture_radius_nm"
            ]
        ),
        "conservative_aperture_diameter_nm": (
            selected[
                "conservative_aperture_diameter_nm"
            ]
        ),
        "open_area_fraction": selected[
            "open_area_fraction"
        ],
        "initial_water_molecules": (
            R0_WATERS
        ),
        "removed_water_molecules": selected[
            "removed_water_molecules"
        ],
        "retained_water_molecules": selected[
            "retained_water_molecules"
        ],
        "initial_lumen_water_molecules": (
            INITIAL_LUMEN_WATERS
        ),
        "removed_lumen_water_molecules": (
            selected[
                "removed_lumen_water_molecules"
            ]
        ),
        "retained_lumen_water_molecules": (
            selected[
                "retained_lumen_water_molecules"
            ]
        ),
        "lumen_water_retention_fraction": (
            selected[
                "lumen_water_retention_fraction"
            ]
        ),
        "minimum_CAP_HBN_distance_nm": (
            selected[
                "minimum_CAP_HBN_distance_nm"
            ]
        ),
        "minimum_CAP_PYR_distance_nm": (
            selected[
                "minimum_CAP_PYR_distance_nm"
            ]
        ),
        "minimum_CAP_OW_distance_nm": (
            selected[
                "minimum_CAP_OW_distance_nm"
            ]
        ),
        "R2_atom_count": (
            expected_selected_atoms
        ),
        "molecule_counts": {
            "HBN": 1,
            "PYR": 4,
            "SOL": selected[
                "retained_water_molecules"
            ],
            "CAPL": 1,
            "CAPU": 1,
        },
        "selected_caps_GRO": relative(
            SELECTED_CAPS_GRO
        ),
        "selected_system_GRO": relative(
            SELECTED_SYSTEM_GRO
        ),
        "topology_generation_authorized": True,
        "energy_minimization_authorized": False,
        "MD_execution_authorized": False,
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "required_next_step": (
            "BUILD_R2_TOPOLOGY_AND_RUN_STATIC_CAP_WATER_SCAN"
        ),
    }

    SELECTED_DEFINITION_JSON.write_text(
        json.dumps(
            selected_definition,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    global_gates = {
        "R1_positive_control_is_validated": (
            r1_summary.get(
                "decision"
            )
            ==
            "R1_FROZEN_SOLUTE_50PS_POSITIVE_CONTROL_VALIDATED"
        ),
        "R1_authorized_R2_static_design": (
            parse_bool(
                r1_summary.get(
                    "R2_static_design_authorized",
                    "false",
                )
            )
        ),
        "R0_atom_count_is_68320": (
            len(r0_atoms)
            == EXPECTED_R0_ATOMS
        ),
        "R1_cap_lattice_has_326_beads": (
            len(cap_atoms)
            == R1_TOTAL_CAP_BEADS
        ),
        "R0_initial_lumen_occupancy_is_428": (
            initial_lumen_count
            == INITIAL_LUMEN_WATERS
        ),
        "at_least_one_candidate_passed": (
            len(
                passing_candidates
            )
            > 0
        ),
        "selected_candidate_passed_all_static_gates": (
            bool(
                selected[
                    "candidate_pass"
                ]
            )
        ),
        "selected_caps_are_symmetric": (
            selected[
                "cap_beads_lower"
            ]
            ==
            selected[
                "cap_beads_upper"
            ]
        ),
        "selected_system_atom_count_is_consistent": (
            len(selected_atoms)
            == expected_selected_atoms
        ),
    }

    failed_global_gates = [
        name
        for name, passed
        in global_gates.items()
        if not passed
    ]

    accepted = (
        len(
            failed_global_gates
        )
        == 0
    )

    decision = (
        "R2_PARTIAL_CAP_GEOMETRY_STATIC_GATE_PASSED"
        if accepted
        else
        "R2_PARTIAL_CAP_GEOMETRY_REQUIRES_REVIEW"
    )

    required_next_step = (
        "BUILD_R2_TOPOLOGY_AND_RUN_STATIC_CAP_WATER_SCAN"
        if accepted
        else
        "REVIEW_R2_GEOMETRY_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "candidate_count": (
            len(candidates)
        ),
        "passing_candidate_count": (
            len(
                passing_candidates
            )
        ),
        "selected_candidate_id": (
            selected[
                "candidate_id"
            ]
        ),
        "selected_removal_radius_nm": (
            selected[
                "removal_radius_nm"
            ]
        ),
        "selected_cap_beads_lower": (
            selected[
                "cap_beads_lower"
            ]
        ),
        "selected_cap_beads_upper": (
            selected[
                "cap_beads_upper"
            ]
        ),
        "selected_cap_beads_total": (
            selected[
                "cap_beads_total"
            ]
        ),
        "selected_removed_cap_beads_total": (
            selected[
                "removed_cap_beads_total"
            ]
        ),
        "selected_effective_aperture_radius_5kBT_nm": (
            selected[
                "effective_aperture_radius_5kBT_nm"
            ]
        ),
        "selected_effective_aperture_diameter_5kBT_nm": (
            selected[
                "effective_aperture_diameter_5kBT_nm"
            ]
        ),
        "selected_conservative_aperture_radius_nm": (
            selected[
                "conservative_aperture_radius_nm"
            ]
        ),
        "selected_conservative_aperture_diameter_nm": (
            selected[
                "conservative_aperture_diameter_nm"
            ]
        ),
        "selected_open_area_fraction": (
            selected[
                "open_area_fraction"
            ]
        ),
        "selected_removed_water_molecules": (
            selected[
                "removed_water_molecules"
            ]
        ),
        "selected_retained_water_molecules": (
            selected[
                "retained_water_molecules"
            ]
        ),
        "selected_retained_lumen_water_molecules": (
            selected[
                "retained_lumen_water_molecules"
            ]
        ),
        "selected_lumen_water_retention_fraction": (
            selected[
                "lumen_water_retention_fraction"
            ]
        ),
        "selected_minimum_CAP_HBN_distance_nm": (
            selected[
                "minimum_CAP_HBN_distance_nm"
            ]
        ),
        "selected_minimum_CAP_PYR_distance_nm": (
            selected[
                "minimum_CAP_PYR_distance_nm"
            ]
        ),
        "selected_minimum_CAP_OW_distance_nm": (
            selected[
                "minimum_CAP_OW_distance_nm"
            ]
        ),
        "selected_R2_atom_count": (
            expected_selected_atoms
        ),
        "failed_gates": (
            " | ".join(
                failed_global_gates
            )
        ),
        "topology_generation_authorized": (
            accepted
        ),
        "energy_minimization_authorized": False,
        "MD_execution_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    write_csv(
        GATE_CSV,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed
            in global_gates.items()
        ],
    )

    candidate_lines = "\n".join(
        (
            f"- `{candidate['candidate_id']}`: "
            f"remove radius="
            f"{candidate['removal_radius_nm']:.3f} nm; "
            f"caps="
            f"{candidate['cap_beads_lower']}/"
            f"{candidate['cap_beads_upper']}; "
            f"effective aperture radius="
            f"{candidate['effective_aperture_radius_5kBT_nm']:.6f} nm; "
            f"open area="
            f"{candidate['open_area_fraction']:.6f}; "
            f"retained lumen waters="
            f"{candidate['retained_lumen_water_molecules']}/428; "
            f"status="
            f"{'PASS' if candidate['candidate_pass'] else 'FAIL'}"
        )
        for candidate in candidates
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in global_gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Partial-Cap Geometry Design

## Scope

R2 introduces symmetric central axial apertures into the validated R1
steric-cap lattice.

R2 remains a neutral frozen steric screening architecture. This gate
does not claim chemical realizability and does not authorize molecular
dynamics.

## Design basis

- Authoritative hydrated source:
  `{relative(R0_T0_GRO)}`
- Validated R1 cap lattice:
  `{relative(R1_CAPS_GRO)}`
- Initial lumen occupancy:
  **{initial_lumen_count} waters**
- R1 cap beads:
  **{R1_CAP_BEADS_PER_END} per end**
- CAP–OW 5 kBT distance:
  **{CAP_WATER_5KBT_DISTANCE_NM:.3f} nm**
- Initial-overlap cutoff:
  **{WATER_OVERLAP_CUTOFF_NM:.3f} nm**
- Target effective aperture radius:
  **{TARGET_EFFECTIVE_APERTURE_RADIUS_NM:.3f} nm**

## Candidate scan

{candidate_lines}

## Selected R2 candidate

- Candidate:
  **{selected['candidate_id']}**
- Nominal cap-bead removal radius:
  **{selected['removal_radius_nm']:.6f} nm**
- Lower/upper cap beads:
  **{selected['cap_beads_lower']}/
  {selected['cap_beads_upper']}**
- Removed cap beads:
  **{selected['removed_cap_beads_total']}**
- Effective 5 kBT aperture radius/diameter:
  **{selected['effective_aperture_radius_5kBT_nm']:.6f}/
  {selected['effective_aperture_diameter_5kBT_nm']:.6f} nm**
- Conservative aperture radius/diameter:
  **{selected['conservative_aperture_radius_nm']:.6f}/
  {selected['conservative_aperture_diameter_nm']:.6f} nm**
- Open-area fraction:
  **{selected['open_area_fraction']:.6f}**
- Removed/retained water molecules:
  **{selected['removed_water_molecules']}/
  {selected['retained_water_molecules']}**
- Retained lumen waters:
  **{selected['retained_lumen_water_molecules']}/
  {INITIAL_LUMEN_WATERS}**
- Lumen-water retention fraction:
  **{selected['lumen_water_retention_fraction']:.6f}**
- Minimum CAP–HBN distance:
  **{selected['minimum_CAP_HBN_distance_nm']:.6f} nm**
- Minimum CAP–PYR distance:
  **{selected['minimum_CAP_PYR_distance_nm']:.6f} nm**
- Minimum CAP–OW distance:
  **{selected['minimum_CAP_OW_distance_nm']:.6f} nm**
- R2 atom count:
  **{expected_selected_atoms}**

## Static gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_global_gates else ' | '.join(failed_global_gates)}**
- Topology generation authorized:
  **{'YES' if accepted else 'NO'}**
- Energy minimization authorized:
  **NO**
- Molecular dynamics authorized:
  **NO**
- Required next step:
  `{required_next_step}`

The selected aperture is a geometry-screening target. Actual water
exchange and retention must be established through a subsequent
validated topology, static interaction scan, water-only minimization,
and short frozen-solute trajectory.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R2 partial-cap geometry design "
        "and static validation completed."
    )

    print(
        "Candidates evaluated / passing: "
        f"{len(candidates)}/"
        f"{len(passing_candidates)}"
    )

    print(
        "Selected candidate / removal radius: "
        f"{selected['candidate_id']} / "
        f"{selected['removal_radius_nm']:.6f} nm"
    )

    print(
        "Selected lower/upper/total cap beads: "
        f"{selected['cap_beads_lower']}/"
        f"{selected['cap_beads_upper']}/"
        f"{selected['cap_beads_total']}"
    )

    print(
        "Removed cap beads: "
        f"{selected['removed_cap_beads_total']}"
    )

    print(
        "Effective aperture radius/diameter at 5 kBT: "
        f"{selected['effective_aperture_radius_5kBT_nm']:.6f}/"
        f"{selected['effective_aperture_diameter_5kBT_nm']:.6f} nm"
    )

    print(
        "Conservative aperture radius/diameter: "
        f"{selected['conservative_aperture_radius_nm']:.6f}/"
        f"{selected['conservative_aperture_diameter_nm']:.6f} nm"
    )

    print(
        "Open-area fraction: "
        f"{selected['open_area_fraction']:.6f}"
    )

    print(
        "Removed / retained waters: "
        f"{selected['removed_water_molecules']}/"
        f"{selected['retained_water_molecules']}"
    )

    print(
        "Retained lumen waters / fraction: "
        f"{selected['retained_lumen_water_molecules']}/"
        f"{INITIAL_LUMEN_WATERS} / "
        f"{selected['lumen_water_retention_fraction']:.6f}"
    )

    print(
        "Minimum CAP-HBN / CAP-PYR / CAP-OW distances: "
        f"{selected['minimum_CAP_HBN_distance_nm']:.6f}/"
        f"{selected['minimum_CAP_PYR_distance_nm']:.6f}/"
        f"{selected['minimum_CAP_OW_distance_nm']:.6f} nm"
    )

    print(
        "Selected R2 atom count: "
        f"{expected_selected_atoms}"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed_global_gates
            else " | ".join(
                failed_global_gates
            )
        )
    )

    print(
        "Topology generation authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Energy minimization authorized: NO"
    )

    print(
        "MD execution authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    print(
        f"Wrote: {relative(CANDIDATE_SCAN_CSV)}"
    )

    print(
        f"Wrote: {relative(SELECTED_CAP_COORDINATES_CSV)}"
    )

    print(
        f"Wrote: {relative(SELECTED_REMOVED_WATERS_CSV)}"
    )

    print(
        f"Wrote: {relative(SELECTED_CAPS_GRO)}"
    )

    print(
        f"Wrote: {relative(SELECTED_SYSTEM_GRO)}"
    )

    print(
        f"Wrote: {relative(SELECTED_DEFINITION_JSON)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(GATE_CSV)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if not accepted:
        raise RuntimeError(
            "R2 partial-cap geometry requires review."
        )


if __name__ == "__main__":
    main()
