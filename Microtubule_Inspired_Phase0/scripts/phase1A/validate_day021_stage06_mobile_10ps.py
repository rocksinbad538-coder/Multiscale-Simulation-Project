#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

PROTOCOL_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

RUN05 = (
    PROTOCOL_ROOT
    / "execution/05_nvt_unrestrained_2ps"
)

RUN06 = (
    PROTOCOL_ROOT
    / "execution/06_nvt_unrestrained_10ps"
)

PREVIOUS_GRO = RUN05 / "05_nvt_unrestrained_2ps.gro"
CURRENT_GRO = RUN06 / "06_nvt_unrestrained_10ps.gro"

OPERATIONAL_SUMMARY = (
    RUN06
    / "06_nvt_unrestrained_10ps_summary.csv"
)

HBN_STRUCTURAL_SUMMARY = (
    RUN06
    / "06_nvt_unrestrained_10ps_structural_summary.csv"
)

HBN_ITP = (
    PROTOCOL_ROOT
    / "protocol_inputs/topology/hbn_bonded_mobile_release.itp"
)

PYR_ITP = (
    PROTOCOL_ROOT
    / "protocol_inputs/topology/pyrene_mobile_release.itp"
)

SUMMARY_CSV = (
    RUN06
    / "stage06_integrated_mobile_pilot_summary.csv"
)

PYR_CSV = (
    RUN06
    / "stage06_pyrene_geometry_summary.csv"
)

REPORT_MD = (
    RUN06
    / "STAGE06_INTEGRATED_MOBILE_PILOT_VALIDATION_DAY021.md"
)

EXPECTED_ATOMS = 68320

HBN_FIRST = 1
HBN_LAST = 1680

PYR_FIRST = 1681
PYR_COUNT = 4
PYR_ATOMS_PER_MOLECULE = 26

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(PROJECT_ROOT)
        )
    except ValueError:
        return str(path.resolve())


def require_inputs() -> None:
    required = (
        PREVIOUS_GRO,
        CURRENT_GRO,
        OPERATIONAL_SUMMARY,
        HBN_STRUCTURAL_SUMMARY,
        HBN_ITP,
        PYR_ITP,
    )

    missing = [
        path
        for path in required
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty required inputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


def read_single_row(
    path: Path,
) -> dict[str, str]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected exactly one row in {path}"
        )

    return rows[0]


def read_gro(
    path: Path,
) -> tuple[list[dict[str, object]], np.ndarray]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    natoms = int(lines[1].strip())

    if natoms != EXPECTED_ATOMS:
        raise RuntimeError(
            f"Unexpected atom count in {path}: {natoms}"
        )

    atom_lines = lines[2 : 2 + natoms]

    atoms: list[dict[str, object]] = []

    for global_index, line in enumerate(
        atom_lines,
        start=1,
    ):
        atoms.append(
            {
                "global_index": global_index,
                "residue_name": line[5:10].strip(),
                "atom_name": line[10:15].strip(),
                "gro_atom_number": int(line[15:20]),
                "position": np.array(
                    [
                        float(line[20:28]),
                        float(line[28:36]),
                        float(line[36:44]),
                    ],
                    dtype=float,
                ),
            }
        )

    box_fields = [
        float(value)
        for value in lines[2 + natoms].split()
    ]

    if len(box_fields) < 3:
        raise RuntimeError(
            f"Invalid simulation box in {path}"
        )

    return (
        atoms,
        np.array(
            box_fields[:3],
            dtype=float,
        ),
    )


def minimum_image(
    vectors: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    return (
        vectors
        - box
        * np.round(
            vectors / box
        )
    )


def parse_itp(
    path: Path,
) -> dict[str, list[tuple]]:
    section = ""

    result: dict[str, list[tuple]] = {
        "bonds": [],
        "angles": [],
        "dihedrals": [],
    }

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line:
            continue

        section_match = SECTION_PATTERN.match(line)

        if section_match:
            section = (
                section_match.group(1)
                .strip()
                .lower()
            )
            continue

        if line.startswith("#"):
            continue

        fields = line.split()

        try:
            if section == "bonds" and len(fields) >= 3:
                equilibrium = (
                    float(fields[3])
                    if len(fields) >= 5
                    else math.nan
                )

                result["bonds"].append(
                    (
                        int(fields[0]),
                        int(fields[1]),
                        int(fields[2]),
                        equilibrium,
                    )
                )

            elif section == "angles" and len(fields) >= 4:
                equilibrium = (
                    float(fields[4])
                    if len(fields) >= 6
                    else math.nan
                )

                result["angles"].append(
                    (
                        int(fields[0]),
                        int(fields[1]),
                        int(fields[2]),
                        int(fields[3]),
                        equilibrium,
                    )
                )

            elif section == "dihedrals" and len(fields) >= 5:
                equilibrium = (
                    float(fields[5])
                    if len(fields) >= 7
                    else math.nan
                )

                result["dihedrals"].append(
                    (
                        int(fields[0]),
                        int(fields[1]),
                        int(fields[2]),
                        int(fields[3]),
                        int(fields[4]),
                        equilibrium,
                    )
                )

        except ValueError:
            continue

    return result


def unwrap_by_bonds(
    wrapped: np.ndarray,
    box: np.ndarray,
    bonds: list[tuple],
) -> np.ndarray:
    count = len(wrapped)

    adjacency: list[list[int]] = [
        []
        for _ in range(count)
    ]

    for bond in bonds:
        atom_i = int(bond[0]) - 1
        atom_j = int(bond[1]) - 1

        if not (
            0 <= atom_i < count
            and 0 <= atom_j < count
        ):
            continue

        adjacency[atom_i].append(atom_j)
        adjacency[atom_j].append(atom_i)

    unwrapped = np.zeros_like(wrapped)
    visited = np.zeros(count, dtype=bool)

    for root in range(count):
        if visited[root]:
            continue

        visited[root] = True
        unwrapped[root] = wrapped[root]
        queue = [root]

        while queue:
            atom_i = queue.pop(0)

            for atom_j in adjacency[atom_i]:
                if visited[atom_j]:
                    continue

                delta = minimum_image(
                    (
                        wrapped[atom_j]
                        - wrapped[atom_i]
                    )[None, :],
                    box,
                )[0]

                unwrapped[atom_j] = (
                    unwrapped[atom_i]
                    + delta
                )

                visited[atom_j] = True
                queue.append(atom_j)

    return unwrapped


def align_current_to_previous(
    previous_wrapped: np.ndarray,
    previous_unwrapped: np.ndarray,
    current_wrapped: np.ndarray,
    current_unwrapped: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    anchor_delta = minimum_image(
        (
            current_wrapped[0]
            - previous_wrapped[0]
        )[None, :],
        box,
    )[0]

    target_anchor = (
        previous_unwrapped[0]
        + anchor_delta
    )

    return (
        current_unwrapped
        + target_anchor
        - current_unwrapped[0]
    )


def kabsch_metrics(
    reference: np.ndarray,
    target: np.ndarray,
) -> tuple[np.ndarray, float]:
    reference_centered = (
        reference
        - reference.mean(axis=0)
    )

    target_centered = (
        target
        - target.mean(axis=0)
    )

    covariance = (
        target_centered.T
        @ reference_centered
    )

    u_matrix, _, vt_matrix = np.linalg.svd(
        covariance
    )

    correction = np.eye(3)

    correction[2, 2] = np.sign(
        np.linalg.det(
            u_matrix @ vt_matrix
        )
    )

    rotation = (
        u_matrix
        @ correction
        @ vt_matrix
    )

    aligned = (
        target_centered
        @ rotation
    )

    residuals = np.linalg.norm(
        aligned - reference_centered,
        axis=1,
    )

    cosine = (
        np.trace(rotation) - 1.0
    ) / 2.0

    cosine = max(
        -1.0,
        min(
            1.0,
            float(cosine),
        ),
    )

    rotation_angle = math.degrees(
        math.acos(cosine)
    )

    return residuals, rotation_angle


def bond_lengths(
    positions: np.ndarray,
    bonds: list[tuple],
) -> np.ndarray:
    return np.array(
        [
            np.linalg.norm(
                positions[int(bond[1]) - 1]
                - positions[int(bond[0]) - 1]
            )
            for bond in bonds
        ],
        dtype=float,
    )


def angle_values(
    positions: np.ndarray,
    angles: list[tuple],
) -> np.ndarray:
    values = []

    for angle in angles:
        atom_i = int(angle[0]) - 1
        atom_j = int(angle[1]) - 1
        atom_k = int(angle[2]) - 1

        vector_ji = (
            positions[atom_i]
            - positions[atom_j]
        )

        vector_jk = (
            positions[atom_k]
            - positions[atom_j]
        )

        denominator = (
            np.linalg.norm(vector_ji)
            * np.linalg.norm(vector_jk)
        )

        if denominator <= 0.0:
            raise RuntimeError(
                "Zero-length vector in angle calculation"
            )

        cosine = float(
            np.dot(vector_ji, vector_jk)
            / denominator
        )

        cosine = max(
            -1.0,
            min(
                1.0,
                cosine,
            ),
        )

        values.append(
            math.degrees(
                math.acos(cosine)
            )
        )

    return np.array(
        values,
        dtype=float,
    )


def dihedral_values(
    positions: np.ndarray,
    box: np.ndarray,
    dihedrals: list[tuple],
) -> np.ndarray:
    values = []

    for dihedral in dihedrals:
        atom_i = int(dihedral[0]) - 1
        atom_j = int(dihedral[1]) - 1
        atom_k = int(dihedral[2]) - 1
        atom_l = int(dihedral[3]) - 1

        b0 = minimum_image(
            (
                positions[atom_j]
                - positions[atom_i]
            )[None, :],
            box,
        )[0]

        b1 = minimum_image(
            (
                positions[atom_k]
                - positions[atom_j]
            )[None, :],
            box,
        )[0]

        b2 = minimum_image(
            (
                positions[atom_l]
                - positions[atom_k]
            )[None, :],
            box,
        )[0]

        b1_norm = np.linalg.norm(b1)

        if b1_norm <= 0.0:
            raise RuntimeError(
                "Zero-length central dihedral bond"
            )

        b1_unit = b1 / b1_norm

        vector_v = (
            b0
            - np.dot(b0, b1_unit)
            * b1_unit
        )

        vector_w = (
            b2
            - np.dot(b2, b1_unit)
            * b1_unit
        )

        x_value = float(
            np.dot(
                vector_v,
                vector_w,
            )
        )

        y_value = float(
            np.dot(
                np.cross(
                    b1_unit,
                    vector_v,
                ),
                vector_w,
            )
        )

        values.append(
            math.degrees(
                math.atan2(
                    y_value,
                    x_value,
                )
            )
        )

    return np.array(
        values,
        dtype=float,
    )


def circular_difference(
    first: np.ndarray,
    second: np.ndarray,
) -> np.ndarray:
    return np.abs(
        (
            first
            - second
            + 180.0
        )
        % 360.0
        - 180.0
    )


def finite_equilibria(
    records: list[tuple],
) -> np.ndarray:
    return np.array(
        [
            float(record[-1])
            for record in records
        ],
        dtype=float,
    )


def pair_contact_metrics(
    positions_a: np.ndarray,
    positions_b: np.ndarray,
    box: np.ndarray,
    thresholds: tuple[float, ...] = (
        0.10,
        0.12,
        0.15,
        0.20,
    ),
    block_size: int = 512,
) -> dict[str, float | int]:
    minimum = math.inf

    counts = {
        threshold: 0
        for threshold in thresholds
    }

    for start in range(
        0,
        len(positions_b),
        block_size,
    ):
        block = positions_b[
            start : start + block_size
        ]

        vectors = (
            positions_a[:, None, :]
            - block[None, :, :]
        )

        vectors = minimum_image(
            vectors,
            box,
        )

        distances = np.linalg.norm(
            vectors,
            axis=2,
        )

        minimum = min(
            minimum,
            float(distances.min()),
        )

        for threshold in thresholds:
            counts[threshold] += int(
                np.count_nonzero(
                    distances < threshold
                )
            )

    result: dict[str, float | int] = {
        "minimum_nm": minimum,
    }

    for threshold in thresholds:
        label = str(threshold).replace(
            ".",
            "p",
        )

        result[
            f"pairs_below_{label}_nm"
        ] = counts[threshold]

    return result


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fieldnames: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for key in row:
            if key in seen:
                continue

            seen.add(key)
            fieldnames.append(key)

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

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fieldnames
                }
            )


def quantile(
    values: np.ndarray,
    probability: float,
) -> float:
    if len(values) == 0:
        return math.nan

    return float(
        np.quantile(
            values,
            probability,
        )
    )


def main() -> None:
    require_inputs()

    operational = read_single_row(
        OPERATIONAL_SUMMARY
    )

    hbn_structural = read_single_row(
        HBN_STRUCTURAL_SUMMARY
    )

    operational_pass = (
        operational.get(
            "decision",
            "",
        ).strip().upper()
        == "PASS"
    )

    hbn_pass = (
        hbn_structural.get(
            "structural_screen",
            "",
        ).strip().upper()
        == "STABLE_CANDIDATE"
    )

    previous_atoms, previous_box = read_gro(
        PREVIOUS_GRO
    )

    current_atoms, current_box = read_gro(
        CURRENT_GRO
    )

    if not np.allclose(
        previous_box,
        current_box,
        atol=1.0e-8,
    ):
        raise RuntimeError(
            "Simulation box changed during Stage06"
        )

    box = current_box

    hbn_topology = parse_itp(
        HBN_ITP
    )

    pyr_topology = parse_itp(
        PYR_ITP
    )

    if not pyr_topology["bonds"]:
        raise RuntimeError(
            "No PYR bonds were parsed"
        )

    if not pyr_topology["angles"]:
        raise RuntimeError(
            "No PYR angles were parsed"
        )

    improper_records = [
        record
        for record in hbn_topology[
            "dihedrals"
        ]
        if int(record[4]) in {
            2,
            4,
        }
    ]

    if not improper_records:
        raise RuntimeError(
            "No HBN improper dihedrals were identified"
        )

    pyr_rows: list[
        dict[str, object]
    ] = []

    all_pyr_bond_deviations = []
    all_pyr_angle_deviations = []

    max_pyr_aligned_rms = 0.0
    max_pyr_aligned_atom = 0.0

    for molecule_index in range(
        PYR_COUNT
    ):
        first_global = (
            PYR_FIRST
            + molecule_index
            * PYR_ATOMS_PER_MOLECULE
        )

        last_global = (
            first_global
            + PYR_ATOMS_PER_MOLECULE
            - 1
        )

        previous_wrapped = np.array(
            [
                previous_atoms[index - 1][
                    "position"
                ]
                for index in range(
                    first_global,
                    last_global + 1,
                )
            ],
            dtype=float,
        )

        current_wrapped = np.array(
            [
                current_atoms[index - 1][
                    "position"
                ]
                for index in range(
                    first_global,
                    last_global + 1,
                )
            ],
            dtype=float,
        )

        previous_unwrapped = unwrap_by_bonds(
            previous_wrapped,
            box,
            pyr_topology["bonds"],
        )

        current_unwrapped = unwrap_by_bonds(
            current_wrapped,
            box,
            pyr_topology["bonds"],
        )

        current_shifted = align_current_to_previous(
            previous_wrapped,
            previous_unwrapped,
            current_wrapped,
            current_unwrapped,
            box,
        )

        residuals, rotation_angle = (
            kabsch_metrics(
                previous_unwrapped,
                current_shifted,
            )
        )

        aligned_rms = float(
            np.sqrt(
                np.mean(
                    residuals ** 2
                )
            )
        )

        aligned_max = float(
            residuals.max()
        )

        max_pyr_aligned_rms = max(
            max_pyr_aligned_rms,
            aligned_rms,
        )

        max_pyr_aligned_atom = max(
            max_pyr_aligned_atom,
            aligned_max,
        )

        centroid_displacement = float(
            np.linalg.norm(
                current_shifted.mean(axis=0)
                - previous_unwrapped.mean(axis=0)
            )
        )

        current_bonds = bond_lengths(
            current_unwrapped,
            pyr_topology["bonds"],
        )

        bond_equilibrium = finite_equilibria(
            pyr_topology["bonds"]
        )

        if np.all(
            np.isfinite(
                bond_equilibrium
            )
        ):
            bond_deviation = np.abs(
                current_bonds
                - bond_equilibrium
            )
        else:
            previous_bonds = bond_lengths(
                previous_unwrapped,
                pyr_topology["bonds"],
            )

            bond_deviation = np.abs(
                current_bonds
                - previous_bonds
            )

        current_angles = angle_values(
            current_unwrapped,
            pyr_topology["angles"],
        )

        angle_equilibrium = finite_equilibria(
            pyr_topology["angles"]
        )

        if np.all(
            np.isfinite(
                angle_equilibrium
            )
        ):
            angle_deviation = np.abs(
                current_angles
                - angle_equilibrium
            )
        else:
            previous_angles = angle_values(
                previous_unwrapped,
                pyr_topology["angles"],
            )

            angle_deviation = np.abs(
                current_angles
                - previous_angles
            )

        all_pyr_bond_deviations.extend(
            bond_deviation.tolist()
        )

        all_pyr_angle_deviations.extend(
            angle_deviation.tolist()
        )

        pyr_rows.append(
            {
                "pyrene": molecule_index + 1,
                "first_global_atom": first_global,
                "last_global_atom": last_global,
                "centroid_displacement_nm": (
                    centroid_displacement
                ),
                "rotation_angle_deg": (
                    rotation_angle
                ),
                "aligned_rms_nm": aligned_rms,
                "aligned_max_nm": aligned_max,
                "bond_deviation_q95_nm": (
                    quantile(
                        bond_deviation,
                        0.95,
                    )
                ),
                "bond_deviation_max_nm": (
                    float(
                        bond_deviation.max()
                    )
                ),
                "angle_deviation_q95_deg": (
                    quantile(
                        angle_deviation,
                        0.95,
                    )
                ),
                "angle_deviation_max_deg": (
                    float(
                        angle_deviation.max()
                    )
                ),
            }
        )

    pyr_bond_deviations = np.array(
        all_pyr_bond_deviations,
        dtype=float,
    )

    pyr_angle_deviations = np.array(
        all_pyr_angle_deviations,
        dtype=float,
    )

    previous_hbn = np.array(
        [
            previous_atoms[index - 1][
                "position"
            ]
            for index in range(
                HBN_FIRST,
                HBN_LAST + 1,
            )
        ],
        dtype=float,
    )

    current_hbn = np.array(
        [
            current_atoms[index - 1][
                "position"
            ]
            for index in range(
                HBN_FIRST,
                HBN_LAST + 1,
            )
        ],
        dtype=float,
    )

    previous_impropers = dihedral_values(
        previous_hbn,
        box,
        improper_records,
    )

    current_impropers = dihedral_values(
        current_hbn,
        box,
        improper_records,
    )

    improper_stage_change = circular_difference(
        current_impropers,
        previous_impropers,
    )

    improper_equilibrium = finite_equilibria(
        improper_records
    )

    if np.all(
        np.isfinite(
            improper_equilibrium
        )
    ):
        improper_equilibrium_deviation = (
            circular_difference(
                current_impropers,
                improper_equilibrium,
            )
        )
    else:
        improper_equilibrium_deviation = (
            improper_stage_change.copy()
        )

    current_positions = np.array(
        [
            atom["position"]
            for atom in current_atoms
        ],
        dtype=float,
    )

    hbn_positions = current_positions[
        HBN_FIRST - 1 : HBN_LAST
    ]

    pyr_indices = np.arange(
        PYR_FIRST - 1,
        PYR_FIRST - 1
        + PYR_COUNT
        * PYR_ATOMS_PER_MOLECULE,
    )

    pyr_heavy_indices = np.array(
        [
            index
            for index in pyr_indices
            if not str(
                current_atoms[index][
                    "atom_name"
                ]
            ).upper().startswith("H")
        ],
        dtype=int,
    )

    water_oxygen_indices = np.array(
        [
            index
            for index, atom in enumerate(
                current_atoms
            )
            if (
                str(
                    atom["residue_name"]
                ).upper()
                == "SOL"
                and str(
                    atom["atom_name"]
                ).upper().startswith("O")
            )
        ],
        dtype=int,
    )

    if len(pyr_heavy_indices) == 0:
        raise RuntimeError(
            "No PYR heavy atoms were identified"
        )

    if len(water_oxygen_indices) == 0:
        raise RuntimeError(
            "No water oxygen atoms were identified"
        )

    hbn_pyr_contacts = pair_contact_metrics(
        hbn_positions,
        current_positions[
            pyr_heavy_indices
        ],
        box,
    )

    pyr_water_contacts = pair_contact_metrics(
        current_positions[
            pyr_heavy_indices
        ],
        current_positions[
            water_oxygen_indices
        ],
        box,
    )

    pyr_pyr_minimum = math.inf
    pyr_pyr_counts = {
        0.10: 0,
        0.12: 0,
        0.15: 0,
        0.20: 0,
    }

    for first_molecule in range(
        PYR_COUNT
    ):
        first_start = (
            PYR_FIRST - 1
            + first_molecule
            * PYR_ATOMS_PER_MOLECULE
        )

        first_stop = (
            first_start
            + PYR_ATOMS_PER_MOLECULE
        )

        first_indices = [
            index
            for index in range(
                first_start,
                first_stop,
            )
            if not str(
                current_atoms[index][
                    "atom_name"
                ]
            ).upper().startswith("H")
        ]

        for second_molecule in range(
            first_molecule + 1,
            PYR_COUNT,
        ):
            second_start = (
                PYR_FIRST - 1
                + second_molecule
                * PYR_ATOMS_PER_MOLECULE
            )

            second_stop = (
                second_start
                + PYR_ATOMS_PER_MOLECULE
            )

            second_indices = [
                index
                for index in range(
                    second_start,
                    second_stop,
                )
                if not str(
                    current_atoms[index][
                        "atom_name"
                    ]
                ).upper().startswith("H")
            ]

            metrics = pair_contact_metrics(
                current_positions[
                    first_indices
                ],
                current_positions[
                    second_indices
                ],
                box,
            )

            pyr_pyr_minimum = min(
                pyr_pyr_minimum,
                float(
                    metrics["minimum_nm"]
                ),
            )

            for threshold in pyr_pyr_counts:
                label = str(
                    threshold
                ).replace(
                    ".",
                    "p",
                )

                pyr_pyr_counts[
                    threshold
                ] += int(
                    metrics[
                        f"pairs_below_{label}_nm"
                    ]
                )

    minimum_intergroup_contact = min(
        float(
            hbn_pyr_contacts[
                "minimum_nm"
            ]
        ),
        pyr_pyr_minimum,
        float(
            pyr_water_contacts[
                "minimum_nm"
            ]
        ),
    )

    blocked_reasons: list[str] = []
    review_reasons: list[str] = []

    if not operational_pass:
        blocked_reasons.append(
            "Stage06 operational decision is not PASS"
        )

    if not hbn_pass:
        blocked_reasons.append(
            "HBN structural screen is not STABLE_CANDIDATE"
        )

    if minimum_intergroup_contact < 0.10:
        blocked_reasons.append(
            "intermolecular heavy-atom contact below 0.10 nm"
        )

    if float(
        pyr_bond_deviations.max()
    ) > 0.03:
        blocked_reasons.append(
            "PYR maximum bond deviation exceeds 0.03 nm"
        )

    if float(
        pyr_angle_deviations.max()
    ) > 25.0:
        blocked_reasons.append(
            "PYR maximum angle deviation exceeds 25 degrees"
        )

    if float(
        improper_stage_change.max()
    ) > 60.0:
        blocked_reasons.append(
            "HBN improper change exceeds 60 degrees"
        )

    if (
        minimum_intergroup_contact
        < 0.14
    ):
        review_reasons.append(
            "intermolecular heavy-atom contact below 0.14 nm"
        )

    if quantile(
        pyr_bond_deviations,
        0.99,
    ) > 0.015:
        review_reasons.append(
            "PYR q99 bond deviation exceeds 0.015 nm"
        )

    if quantile(
        pyr_angle_deviations,
        0.99,
    ) > 10.0:
        review_reasons.append(
            "PYR q99 angle deviation exceeds 10 degrees"
        )

    if max_pyr_aligned_rms > 0.08:
        review_reasons.append(
            "a PYR aligned RMS exceeds 0.08 nm"
        )

    if quantile(
        improper_stage_change,
        0.99,
    ) > 20.0:
        review_reasons.append(
            "HBN improper q99 change exceeds 20 degrees"
        )

    if blocked_reasons:
        readiness = "BLOCKED"
        next_step = "REVIEW_MODEL"
    elif review_reasons:
        readiness = "REVIEW"
        next_step = "ADDITIONAL_2PS_OR_TARGETED_DIAGNOSTIC"
    else:
        readiness = (
            "READY_FOR_EXTENDED_MOBILE_VALIDATION"
        )
        next_step = (
            "EXTENDED_MOBILE_VALIDATION_CANDIDATE"
        )

    summary = {
        "stage": "06_nvt_unrestrained_10ps",
        "operational_pass": operational_pass,
        "hbn_structural_screen": (
            hbn_structural.get(
                "structural_screen",
                "",
            )
        ),
        "pyrene_count": PYR_COUNT,
        "PYR_max_aligned_rms_nm": (
            max_pyr_aligned_rms
        ),
        "PYR_max_aligned_atom_residual_nm": (
            max_pyr_aligned_atom
        ),
        "PYR_bond_deviation_q95_nm": (
            quantile(
                pyr_bond_deviations,
                0.95,
            )
        ),
        "PYR_bond_deviation_q99_nm": (
            quantile(
                pyr_bond_deviations,
                0.99,
            )
        ),
        "PYR_bond_deviation_max_nm": (
            float(
                pyr_bond_deviations.max()
            )
        ),
        "PYR_angle_deviation_q95_deg": (
            quantile(
                pyr_angle_deviations,
                0.95,
            )
        ),
        "PYR_angle_deviation_q99_deg": (
            quantile(
                pyr_angle_deviations,
                0.99,
            )
        ),
        "PYR_angle_deviation_max_deg": (
            float(
                pyr_angle_deviations.max()
            )
        ),
        "HBN_improper_count": (
            len(improper_records)
        ),
        "HBN_improper_stage_change_q95_deg": (
            quantile(
                improper_stage_change,
                0.95,
            )
        ),
        "HBN_improper_stage_change_q99_deg": (
            quantile(
                improper_stage_change,
                0.99,
            )
        ),
        "HBN_improper_stage_change_max_deg": (
            float(
                improper_stage_change.max()
            )
        ),
        "HBN_improper_equilibrium_deviation_q99_deg": (
            quantile(
                improper_equilibrium_deviation,
                0.99,
            )
        ),
        "HBN_improper_equilibrium_deviation_max_deg": (
            float(
                improper_equilibrium_deviation.max()
            )
        ),
        "HBN_PYR_minimum_heavy_contact_nm": (
            hbn_pyr_contacts[
                "minimum_nm"
            ]
        ),
        "PYR_PYR_minimum_heavy_contact_nm": (
            pyr_pyr_minimum
        ),
        "PYR_waterO_minimum_contact_nm": (
            pyr_water_contacts[
                "minimum_nm"
            ]
        ),
        "minimum_intergroup_contact_nm": (
            minimum_intergroup_contact
        ),
        "HBN_PYR_pairs_below_0p15_nm": (
            hbn_pyr_contacts[
                "pairs_below_0p15_nm"
            ]
        ),
        "PYR_PYR_pairs_below_0p15_nm": (
            pyr_pyr_counts[0.15]
        ),
        "PYR_waterO_pairs_below_0p15_nm": (
            pyr_water_contacts[
                "pairs_below_0p15_nm"
            ]
        ),
        "pilot_readiness": readiness,
        "authorized_next_step": (
            next_step
        ),
        "long_mobile_production_authorized": (
            False
        ),
        "review_reasons": (
            " | ".join(
                review_reasons
            )
        ),
        "blocked_reasons": (
            " | ".join(
                blocked_reasons
            )
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    write_csv(
        PYR_CSV,
        pyr_rows,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Stage06 Integrated Mobile-Pilot Validation\n\n"
        )

        handle.write(
            f"- Pilot readiness: **{readiness}**\n"
        )

        handle.write(
            f"- Authorized next step: `{next_step}`\n"
        )

        handle.write(
            "- Long mobile production authorized: **NO**\n\n"
        )

        handle.write(
            "## Pyrene geometry\n\n"
        )

        handle.write(
            f"- Maximum aligned RMS across PYR molecules: "
            f"{max_pyr_aligned_rms:.8f} nm\n"
        )

        handle.write(
            f"- Bond deviation q99/max: "
            f"{quantile(pyr_bond_deviations, 0.99):.8f}/"
            f"{pyr_bond_deviations.max():.8f} nm\n"
        )

        handle.write(
            f"- Angle deviation q99/max: "
            f"{quantile(pyr_angle_deviations, 0.99):.4f}/"
            f"{pyr_angle_deviations.max():.4f} degrees\n\n"
        )

        handle.write(
            "## HBN impropers\n\n"
        )

        handle.write(
            f"- Improper count: {len(improper_records)}\n"
        )

        handle.write(
            f"- Stage change q99/max: "
            f"{quantile(improper_stage_change, 0.99):.4f}/"
            f"{improper_stage_change.max():.4f} degrees\n\n"
        )

        handle.write(
            "## Intermolecular contacts\n\n"
        )

        handle.write(
            f"- HBN–PYR minimum heavy-atom distance: "
            f"{hbn_pyr_contacts['minimum_nm']:.8f} nm\n"
        )

        handle.write(
            f"- PYR–PYR minimum heavy-atom distance: "
            f"{pyr_pyr_minimum:.8f} nm\n"
        )

        handle.write(
            f"- PYR–water oxygen minimum distance: "
            f"{pyr_water_contacts['minimum_nm']:.8f} nm\n"
        )

        if review_reasons:
            handle.write(
                "\n## Review reasons\n\n"
            )

            for reason in review_reasons:
                handle.write(
                    f"- {reason}\n"
                )

        if blocked_reasons:
            handle.write(
                "\n## Blocking reasons\n\n"
            )

            for reason in blocked_reasons:
                handle.write(
                    f"- {reason}\n"
                )

    print(
        "Day021 Stage06 integrated mobile-pilot validation completed."
    )

    print(
        f"Operational decision: "
        f"{'PASS' if operational_pass else 'FAIL'}"
    )

    print(
        f"HBN structural screen: "
        f"{hbn_structural.get('structural_screen', 'UNKNOWN')}"
    )

    print(
        f"PYR molecules analyzed: {PYR_COUNT}"
    )

    print(
        "PYR maximum aligned RMS / atom residual: "
        f"{max_pyr_aligned_rms:.8f}/"
        f"{max_pyr_aligned_atom:.8f} nm"
    )

    print(
        "PYR bond deviation q95/q99/max: "
        f"{quantile(pyr_bond_deviations, 0.95):.8f}/"
        f"{quantile(pyr_bond_deviations, 0.99):.8f}/"
        f"{pyr_bond_deviations.max():.8f} nm"
    )

    print(
        "PYR angle deviation q95/q99/max: "
        f"{quantile(pyr_angle_deviations, 0.95):.4f}/"
        f"{quantile(pyr_angle_deviations, 0.99):.4f}/"
        f"{pyr_angle_deviations.max():.4f} deg"
    )

    print(
        f"HBN improper count: "
        f"{len(improper_records)}"
    )

    print(
        "HBN improper stage-change q95/q99/max: "
        f"{quantile(improper_stage_change, 0.95):.4f}/"
        f"{quantile(improper_stage_change, 0.99):.4f}/"
        f"{improper_stage_change.max():.4f} deg"
    )

    print(
        "Minimum contacts HBN-PYR / PYR-PYR / PYR-waterO: "
        f"{hbn_pyr_contacts['minimum_nm']:.8f}/"
        f"{pyr_pyr_minimum:.8f}/"
        f"{pyr_water_contacts['minimum_nm']:.8f} nm"
    )

    print(
        f"Pilot readiness: {readiness}"
    )

    print(
        f"Authorized next step: {next_step}"
    )

    print(
        "Long mobile production authorized: NO"
    )

    if review_reasons:
        print(
            "Review reasons: "
            + " | ".join(
                review_reasons
            )
        )

    if blocked_reasons:
        print(
            "Blocking reasons: "
            + " | ".join(
                blocked_reasons
            )
        )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )


if __name__ == "__main__":
    main()
