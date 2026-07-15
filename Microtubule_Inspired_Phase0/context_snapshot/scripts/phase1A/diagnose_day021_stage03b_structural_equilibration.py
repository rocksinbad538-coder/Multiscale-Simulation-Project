#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

PROTOCOL_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

RUN02 = (
    PROTOCOL_ROOT
    / "execution/02_nvt_k10000_1ps"
)

RUN03 = (
    PROTOCOL_ROOT
    / "execution/03_nvt_k1000_2ps"
)

RUN03B = (
    PROTOCOL_ROOT
    / "execution/03b_nvt_k1000_hold_2ps"
)

ACCEPTED_GRO = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute/"
    "nvt_100ps_frozenSolute.gro"
)

STAGE02_GRO = (
    RUN02
    / "02_nvt_k10000_1ps.gro"
)

STAGE03_GRO = (
    RUN03
    / "03_nvt_k1000_2ps.gro"
)

STAGE03B_GRO = (
    RUN03B
    / "03b_nvt_k1000_hold_2ps.gro"
)

STAGE03B_LOG = (
    RUN03B
    / "03b_nvt_k1000_hold_2ps.log"
)

STAGE03B_CONSOLE = (
    RUN03B
    / "03b_nvt_k1000_hold_2ps_mdrun_console.log"
)

HBN_ITP = (
    PROTOCOL_ROOT
    / "protocol_inputs/topology/"
    "hbn_bonded_mobile_release.itp"
)

SUMMARY_CSV = (
    RUN03B
    / "stage03b_structural_equilibration_summary.csv"
)

OUTLIERS_CSV = (
    RUN03B
    / "stage03b_hbn_displacement_outliers.csv"
)

BOND_CSV = (
    RUN03B
    / "stage03b_hbn_largest_bond_deviations.csv"
)

ANGLE_CSV = (
    RUN03B
    / "stage03b_hbn_largest_angle_deviations.csv"
)

REPORT_MD = (
    RUN03B
    / "STAGE03B_STRUCTURAL_EQUILIBRATION_DAY021.md"
)

HBN_FIRST = 1
HBN_LAST = 1680

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+)\s*\]\s*$"
)

NONFINITE_PATTERN = re.compile(
    r"(?<![A-Za-z])"
    r"(?:nan|[-+]?inf(?:inity)?)"
    r"(?![A-Za-z])",
    re.IGNORECASE,
)

HARMLESS_EPSILON_RF = re.compile(
    r"^\s*epsilon-rf\s*=\s*"
    r"(?:inf|infinity)\s*$",
    re.IGNORECASE,
)

SERIOUS_PATTERN = re.compile(
    r"LINCS WARNING|"
    r"Fatal error|"
    r"SETTLE.*(?:error|failed|cannot)|"
    r"SHAKE.*(?:failed|did not converge)|"
    r"constraint.*(?:error|failed)",
    re.IGNORECASE,
)


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(
                PROJECT_ROOT
            )
        )
    except ValueError:
        return str(path.resolve())


def require_inputs() -> None:
    required = (
        ACCEPTED_GRO,
        STAGE02_GRO,
        STAGE03_GRO,
        STAGE03B_GRO,
        STAGE03B_LOG,
        STAGE03B_CONSOLE,
        HBN_ITP,
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


def read_gro(
    path: Path,
) -> tuple[list[dict[str, object]], np.ndarray]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    natoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + natoms]

    if len(atom_lines) != natoms:
        raise RuntimeError(
            f"GRO atom-count mismatch: {path}"
        )

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
            f"Invalid box in {path}"
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


def hbn_positions(
    atoms: list[dict[str, object]],
) -> np.ndarray:
    return np.array(
        [
            atoms[index - 1]["position"]
            for index in range(
                HBN_FIRST,
                HBN_LAST + 1,
            )
        ],
        dtype=float,
    )


def kabsch_residuals(
    reference: np.ndarray,
    target: np.ndarray,
) -> np.ndarray:
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

    u_matrix, _, vt_matrix = (
        np.linalg.svd(covariance)
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

    return np.linalg.norm(
        aligned
        - reference_centered,
        axis=1,
    )


def displacement_analysis(
    source: np.ndarray,
    target: np.ndarray,
    box: np.ndarray,
) -> dict[str, object]:
    vectors = minimum_image(
        target - source,
        box,
    )

    magnitudes = np.linalg.norm(
        vectors,
        axis=1,
    )

    mean_vector = vectors.mean(
        axis=0
    )

    translation_removed = np.linalg.norm(
        vectors - mean_vector,
        axis=1,
    )

    target_unwrapped = (
        source + vectors
    )

    aligned = kabsch_residuals(
        source,
        target_unwrapped,
    )

    return {
        "vectors": vectors,
        "magnitudes": magnitudes,
        "aligned": aligned,
        "rms": float(
            np.sqrt(
                np.mean(
                    magnitudes ** 2
                )
            )
        ),
        "median": float(
            np.median(magnitudes)
        ),
        "q95": float(
            np.quantile(
                magnitudes,
                0.95,
            )
        ),
        "q99": float(
            np.quantile(
                magnitudes,
                0.99,
            )
        ),
        "maximum": float(
            magnitudes.max()
        ),
        "maximum_local_index": int(
            magnitudes.argmax()
        ),
        "translation": float(
            np.linalg.norm(
                mean_vector
            )
        ),
        "translation_removed_rms": float(
            np.sqrt(
                np.mean(
                    translation_removed ** 2
                )
            )
        ),
        "translation_removed_max": float(
            translation_removed.max()
        ),
        "aligned_rms": float(
            np.sqrt(
                np.mean(
                    aligned ** 2
                )
            )
        ),
        "aligned_max": float(
            aligned.max()
        ),
        "above_0p05": int(
            np.count_nonzero(
                magnitudes > 0.05
            )
        ),
        "above_0p08": int(
            np.count_nonzero(
                magnitudes > 0.08
            )
        ),
        "above_0p10": int(
            np.count_nonzero(
                magnitudes > 0.10
            )
        ),
        "above_0p12": int(
            np.count_nonzero(
                magnitudes > 0.12
            )
        ),
        "above_0p15": int(
            np.count_nonzero(
                magnitudes > 0.15
            )
        ),
    }


def parse_hbn_itp(
    path: Path,
) -> tuple[
    list[tuple[int, int, float]],
    list[tuple[int, int, int, float]],
]:
    current_section = ""

    bonds: list[
        tuple[int, int, float]
    ] = []

    angles: list[
        tuple[int, int, int, float]
    ] = []

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

        section_match = (
            SECTION_PATTERN.match(line)
        )

        if section_match:
            current_section = (
                section_match.group(1)
                .strip()
                .lower()
            )
            continue

        if line.startswith("#"):
            continue

        fields = line.split()

        if current_section == "bonds":
            if len(fields) < 5:
                continue

            try:
                atom_i = int(fields[0])
                atom_j = int(fields[1])
                equilibrium_nm = float(
                    fields[3]
                )
            except ValueError:
                continue

            bonds.append(
                (
                    atom_i,
                    atom_j,
                    equilibrium_nm,
                )
            )

        elif current_section == "angles":
            if len(fields) < 6:
                continue

            try:
                atom_i = int(fields[0])
                atom_j = int(fields[1])
                atom_k = int(fields[2])
                equilibrium_deg = float(
                    fields[4]
                )
            except ValueError:
                continue

            angles.append(
                (
                    atom_i,
                    atom_j,
                    atom_k,
                    equilibrium_deg,
                )
            )

    if not bonds:
        raise RuntimeError(
            "No explicit HBN bonds were parsed"
        )

    if not angles:
        raise RuntimeError(
            "No explicit HBN angles were parsed"
        )

    return bonds, angles


def bond_lengths(
    positions: np.ndarray,
    box: np.ndarray,
    bonds: list[
        tuple[int, int, float]
    ],
) -> np.ndarray:
    vectors = np.array(
        [
            positions[atom_j - 1]
            - positions[atom_i - 1]
            for (
                atom_i,
                atom_j,
                _,
            ) in bonds
        ],
        dtype=float,
    )

    vectors = minimum_image(
        vectors,
        box,
    )

    return np.linalg.norm(
        vectors,
        axis=1,
    )


def angle_values(
    positions: np.ndarray,
    box: np.ndarray,
    angles: list[
        tuple[int, int, int, float]
    ],
) -> np.ndarray:
    values = []

    for (
        atom_i,
        atom_j,
        atom_k,
        _,
    ) in angles:
        vector_ji = minimum_image(
            (
                positions[atom_i - 1]
                - positions[atom_j - 1]
            )[None, :],
            box,
        )[0]

        vector_jk = minimum_image(
            (
                positions[atom_k - 1]
                - positions[atom_j - 1]
            )[None, :],
            box,
        )[0]

        denominator = (
            np.linalg.norm(vector_ji)
            * np.linalg.norm(vector_jk)
        )

        if denominator <= 0.0:
            raise RuntimeError(
                "Zero-length vector in angle evaluation"
            )

        cosine = float(
            np.dot(
                vector_ji,
                vector_jk,
            )
            / denominator
        )

        cosine = min(
            1.0,
            max(
                -1.0,
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


def audit_logs() -> tuple[int, int, int]:
    harmless = 0
    harmful = 0
    serious = 0

    for path in (
        STAGE03B_LOG,
        STAGE03B_CONSOLE,
    ):
        for line in path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines():
            if SERIOUS_PATTERN.search(line):
                serious += 1

            if not NONFINITE_PATTERN.search(
                line
            ):
                continue

            if HARMLESS_EPSILON_RF.match(
                line
            ):
                harmless += 1
            else:
                harmful += 1

    return harmless, harmful, serious


def audit_xvg() -> tuple[int, int]:
    numeric_rows = 0
    nonfinite_values = 0

    for path in sorted(
        RUN03B.glob("*.xvg")
    ):
        for raw_line in path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines():
            line = raw_line.strip()

            if (
                not line
                or line.startswith("#")
                or line.startswith("@")
            ):
                continue

            numeric_rows += 1

            for token in line.split():
                try:
                    value = float(token)
                except ValueError:
                    continue

                if not math.isfinite(value):
                    nonfinite_values += 1

    return numeric_rows, nonfinite_values


def quantile(
    values: np.ndarray,
    probability: float,
) -> float:
    return float(
        np.quantile(
            values,
            probability,
        )
    )


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
                    key: row.get(key, "")
                    for key in fieldnames
                }
            )


def main() -> None:
    require_inputs()

    accepted_atoms, accepted_box = (
        read_gro(
            ACCEPTED_GRO
        )
    )

    stage02_atoms, stage02_box = (
        read_gro(
            STAGE02_GRO
        )
    )

    stage03_atoms, stage03_box = (
        read_gro(
            STAGE03_GRO
        )
    )

    stage03b_atoms, stage03b_box = (
        read_gro(
            STAGE03B_GRO
        )
    )

    atom_counts = {
        len(accepted_atoms),
        len(stage02_atoms),
        len(stage03_atoms),
        len(stage03b_atoms),
    }

    if atom_counts != {68320}:
        raise RuntimeError(
            f"Unexpected atom counts: {atom_counts}"
        )

    for box in (
        stage02_box,
        stage03_box,
        stage03b_box,
    ):
        if not np.allclose(
            box,
            accepted_box,
            atol=1.0e-8,
        ):
            raise RuntimeError(
                "Simulation box changed during NVT"
            )

    accepted_hbn = hbn_positions(
        accepted_atoms
    )

    stage02_hbn = hbn_positions(
        stage02_atoms
    )

    stage03_hbn = hbn_positions(
        stage03_atoms
    )

    stage03b_hbn = hbn_positions(
        stage03b_atoms
    )

    stage03_incremental = (
        displacement_analysis(
            stage02_hbn,
            stage03_hbn,
            stage03_box,
        )
    )

    stage03b_incremental = (
        displacement_analysis(
            stage03_hbn,
            stage03b_hbn,
            stage03b_box,
        )
    )

    stage03b_cumulative = (
        displacement_analysis(
            accepted_hbn,
            stage03b_hbn,
            stage03b_box,
        )
    )

    bonds, angles = parse_hbn_itp(
        HBN_ITP
    )

    bond_equilibrium = np.array(
        [
            equilibrium
            for _, _, equilibrium
            in bonds
        ],
        dtype=float,
    )

    angle_equilibrium = np.array(
        [
            equilibrium
            for _, _, _, equilibrium
            in angles
        ],
        dtype=float,
    )

    stage03_bonds = bond_lengths(
        stage03_hbn,
        stage03_box,
        bonds,
    )

    stage03b_bonds = bond_lengths(
        stage03b_hbn,
        stage03b_box,
        bonds,
    )

    stage03_angles = angle_values(
        stage03_hbn,
        stage03_box,
        angles,
    )

    stage03b_angles = angle_values(
        stage03b_hbn,
        stage03b_box,
        angles,
    )

    stage03b_bond_eq_deviation = np.abs(
        stage03b_bonds
        - bond_equilibrium
    )

    stage03_to_stage03b_bond_change = np.abs(
        stage03b_bonds
        - stage03_bonds
    )

    stage03b_angle_eq_deviation = np.abs(
        stage03b_angles
        - angle_equilibrium
    )

    stage03_to_stage03b_angle_change = np.abs(
        stage03b_angles
        - stage03_angles
    )

    top03 = set(
        np.argsort(
            stage03_incremental[
                "magnitudes"
            ]
        )[::-1][:20].tolist()
    )

    top03b_order = np.argsort(
        stage03b_incremental[
            "magnitudes"
        ]
    )[::-1]

    top03b = set(
        top03b_order[:20].tolist()
    )

    top20_overlap = len(
        top03 & top03b
    )

    outlier_rows = []

    for rank, local_index in enumerate(
        top03b_order[:30],
        start=1,
    ):
        global_index = (
            HBN_FIRST
            + int(local_index)
        )

        atom = stage03b_atoms[
            global_index - 1
        ]

        outlier_rows.append(
            {
                "rank": rank,
                "global_atom_index": (
                    global_index
                ),
                "atom_name": (
                    atom["atom_name"]
                ),
                "residue_name": (
                    atom["residue_name"]
                ),
                "stage03_incremental_nm": float(
                    stage03_incremental[
                        "magnitudes"
                    ][local_index]
                ),
                "stage03b_incremental_nm": float(
                    stage03b_incremental[
                        "magnitudes"
                    ][local_index]
                ),
                "stage03b_cumulative_nm": float(
                    stage03b_cumulative[
                        "magnitudes"
                    ][local_index]
                ),
                "stage03b_aligned_residual_nm": float(
                    stage03b_incremental[
                        "aligned"
                    ][local_index]
                ),
                "was_stage03_top20": (
                    int(local_index)
                    in top03
                ),
            }
        )

    bond_order = np.argsort(
        stage03b_bond_eq_deviation
    )[::-1]

    bond_rows = []

    for rank, bond_index in enumerate(
        bond_order[:30],
        start=1,
    ):
        atom_i, atom_j, equilibrium = (
            bonds[int(bond_index)]
        )

        bond_rows.append(
            {
                "rank": rank,
                "atom_i": atom_i,
                "atom_j": atom_j,
                "equilibrium_nm": (
                    equilibrium
                ),
                "stage03_length_nm": float(
                    stage03_bonds[
                        bond_index
                    ]
                ),
                "stage03b_length_nm": float(
                    stage03b_bonds[
                        bond_index
                    ]
                ),
                "stage03b_equilibrium_deviation_nm": float(
                    stage03b_bond_eq_deviation[
                        bond_index
                    ]
                ),
                "stage03_to_stage03b_change_nm": float(
                    stage03_to_stage03b_bond_change[
                        bond_index
                    ]
                ),
            }
        )

    angle_order = np.argsort(
        stage03b_angle_eq_deviation
    )[::-1]

    angle_rows = []

    for rank, angle_index in enumerate(
        angle_order[:30],
        start=1,
    ):
        (
            atom_i,
            atom_j,
            atom_k,
            equilibrium,
        ) = angles[int(angle_index)]

        angle_rows.append(
            {
                "rank": rank,
                "atom_i": atom_i,
                "atom_j": atom_j,
                "atom_k": atom_k,
                "equilibrium_deg": (
                    equilibrium
                ),
                "stage03_deg": float(
                    stage03_angles[
                        angle_index
                    ]
                ),
                "stage03b_deg": float(
                    stage03b_angles[
                        angle_index
                    ]
                ),
                "stage03b_equilibrium_deviation_deg": float(
                    stage03b_angle_eq_deviation[
                        angle_index
                    ]
                ),
                "stage03_to_stage03b_change_deg": float(
                    stage03_to_stage03b_angle_change[
                        angle_index
                    ]
                ),
            }
        )

    displaced_names = Counter(
        stage03b_atoms[
            HBN_FIRST - 1 + index
        ]["atom_name"]
        for index, value in enumerate(
            stage03b_incremental[
                "magnitudes"
            ]
        )
        if value > 0.08
    )

    harmless, harmful, serious = (
        audit_logs()
    )

    xvg_rows, xvg_nonfinite = (
        audit_xvg()
    )

    blocked_reasons = []
    review_flags = []

    if harmful:
        blocked_reasons.append(
            "harmful non-finite log values"
        )

    if serious:
        blocked_reasons.append(
            "serious instability signatures"
        )

    if xvg_nonfinite:
        blocked_reasons.append(
            "non-finite XVG values"
        )

    if float(
        stage03b_bonds.max()
    ) > 0.20:
        blocked_reasons.append(
            "HBN bond longer than 0.20 nm"
        )

    if (
        stage03b_incremental[
            "aligned_rms"
        ]
        > 0.05
    ):
        review_flags.append(
            "aligned HBN RMS exceeds 0.05 nm"
        )

    if quantile(
        stage03b_bond_eq_deviation,
        0.99,
    ) > 0.015:
        review_flags.append(
            "q99 bond-equilibrium deviation "
            "exceeds 0.015 nm"
        )

    if float(
        stage03b_bond_eq_deviation.max()
    ) > 0.03:
        review_flags.append(
            "maximum bond-equilibrium deviation "
            "exceeds 0.03 nm"
        )

    if quantile(
        stage03b_angle_eq_deviation,
        0.95,
    ) > 20.0:
        review_flags.append(
            "q95 angle-equilibrium deviation "
            "exceeds 20 degrees"
        )

    if float(
        stage03b_angle_eq_deviation.max()
    ) > 60.0:
        review_flags.append(
            "maximum angle-equilibrium deviation "
            "exceeds 60 degrees"
        )

    if blocked_reasons:
        structural_screen = "BLOCKED"
    elif review_flags:
        structural_screen = "REVIEW"
    else:
        structural_screen = (
            "STABLE_CANDIDATE"
        )

    summary = {
        "stage": (
            "03b_nvt_k1000_hold_2ps"
        ),
        "structural_screen": (
            structural_screen
        ),
        "harmless_epsilon_rf_matches": (
            harmless
        ),
        "harmful_nonfinite_matches": (
            harmful
        ),
        "serious_instability_matches": (
            serious
        ),
        "numeric_xvg_rows": (
            xvg_rows
        ),
        "numeric_xvg_nonfinite_values": (
            xvg_nonfinite
        ),
        "HBN_incremental_rms_nm": (
            stage03b_incremental["rms"]
        ),
        "HBN_incremental_median_nm": (
            stage03b_incremental["median"]
        ),
        "HBN_incremental_q95_nm": (
            stage03b_incremental["q95"]
        ),
        "HBN_incremental_q99_nm": (
            stage03b_incremental["q99"]
        ),
        "HBN_incremental_max_nm": (
            stage03b_incremental[
                "maximum"
            ]
        ),
        "HBN_incremental_aligned_rms_nm": (
            stage03b_incremental[
                "aligned_rms"
            ]
        ),
        "HBN_incremental_aligned_max_nm": (
            stage03b_incremental[
                "aligned_max"
            ]
        ),
        "HBN_incremental_translation_nm": (
            stage03b_incremental[
                "translation"
            ]
        ),
        "HBN_atoms_above_0p05_nm": (
            stage03b_incremental[
                "above_0p05"
            ]
        ),
        "HBN_atoms_above_0p08_nm": (
            stage03b_incremental[
                "above_0p08"
            ]
        ),
        "HBN_atoms_above_0p10_nm": (
            stage03b_incremental[
                "above_0p10"
            ]
        ),
        "HBN_atoms_above_0p12_nm": (
            stage03b_incremental[
                "above_0p12"
            ]
        ),
        "HBN_atoms_above_0p15_nm": (
            stage03b_incremental[
                "above_0p15"
            ]
        ),
        "HBN_stage03_stage03b_top20_overlap": (
            top20_overlap
        ),
        "HBN_outlier_atom_names_above_0p08_nm": (
            ";".join(
                f"{name}:{count}"
                for name, count
                in sorted(
                    displaced_names.items()
                )
            )
        ),
        "HBN_bond_count": (
            len(bonds)
        ),
        "HBN_bond_length_mean_nm": float(
            stage03b_bonds.mean()
        ),
        "HBN_bond_length_min_nm": float(
            stage03b_bonds.min()
        ),
        "HBN_bond_length_max_nm": float(
            stage03b_bonds.max()
        ),
        "HBN_bond_equilibrium_deviation_q95_nm": (
            quantile(
                stage03b_bond_eq_deviation,
                0.95,
            )
        ),
        "HBN_bond_equilibrium_deviation_q99_nm": (
            quantile(
                stage03b_bond_eq_deviation,
                0.99,
            )
        ),
        "HBN_bond_equilibrium_deviation_max_nm": float(
            stage03b_bond_eq_deviation.max()
        ),
        "HBN_stage03_to_stage03b_bond_change_q95_nm": (
            quantile(
                stage03_to_stage03b_bond_change,
                0.95,
            )
        ),
        "HBN_stage03_to_stage03b_bond_change_max_nm": float(
            stage03_to_stage03b_bond_change.max()
        ),
        "HBN_angle_count": (
            len(angles)
        ),
        "HBN_angle_equilibrium_deviation_q95_deg": (
            quantile(
                stage03b_angle_eq_deviation,
                0.95,
            )
        ),
        "HBN_angle_equilibrium_deviation_q99_deg": (
            quantile(
                stage03b_angle_eq_deviation,
                0.99,
            )
        ),
        "HBN_angle_equilibrium_deviation_max_deg": float(
            stage03b_angle_eq_deviation.max()
        ),
        "HBN_stage03_to_stage03b_angle_change_q95_deg": (
            quantile(
                stage03_to_stage03b_angle_change,
                0.95,
            )
        ),
        "HBN_stage03_to_stage03b_angle_change_max_deg": float(
            stage03_to_stage03b_angle_change.max()
        ),
        "review_flags": (
            " | ".join(
                review_flags
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
        OUTLIERS_CSV,
        outlier_rows,
    )

    write_csv(
        BOND_CSV,
        bond_rows,
    )

    write_csv(
        ANGLE_CSV,
        angle_rows,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Stage03b Structural Equilibration\n\n"
        )

        handle.write(
            f"- Structural screen: "
            f"**{structural_screen}**\n"
        )

        handle.write(
            f"- Serious instability signatures: "
            f"{serious}\n"
        )

        handle.write(
            f"- Harmful non-finite values: "
            f"{harmful}\n"
        )

        handle.write(
            f"- Non-finite XVG values: "
            f"{xvg_nonfinite}\n\n"
        )

        handle.write(
            "## HBN displacement distribution\n\n"
        )

        handle.write(
            f"- RMS/median/q95/q99/max: "
            f"{stage03b_incremental['rms']:.8f}/"
            f"{stage03b_incremental['median']:.8f}/"
            f"{stage03b_incremental['q95']:.8f}/"
            f"{stage03b_incremental['q99']:.8f}/"
            f"{stage03b_incremental['maximum']:.8f} nm\n"
        )

        handle.write(
            f"- Kabsch-aligned RMS/max: "
            f"{stage03b_incremental['aligned_rms']:.8f}/"
            f"{stage03b_incremental['aligned_max']:.8f} nm\n"
        )

        handle.write(
            f"- Atoms above 0.05/0.08/0.10/0.12/0.15 nm: "
            f"{stage03b_incremental['above_0p05']}/"
            f"{stage03b_incremental['above_0p08']}/"
            f"{stage03b_incremental['above_0p10']}/"
            f"{stage03b_incremental['above_0p12']}/"
            f"{stage03b_incremental['above_0p15']}\n"
        )

        handle.write(
            f"- Stage03/Stage03b top-20 overlap: "
            f"{top20_overlap}/20\n\n"
        )

        handle.write(
            "## HBN bonded geometry\n\n"
        )

        handle.write(
            f"- Bond length mean/min/max: "
            f"{stage03b_bonds.mean():.8f}/"
            f"{stage03b_bonds.min():.8f}/"
            f"{stage03b_bonds.max():.8f} nm\n"
        )

        handle.write(
            f"- Bond equilibrium deviation q95/q99/max: "
            f"{quantile(stage03b_bond_eq_deviation, 0.95):.8f}/"
            f"{quantile(stage03b_bond_eq_deviation, 0.99):.8f}/"
            f"{stage03b_bond_eq_deviation.max():.8f} nm\n"
        )

        handle.write(
            f"- Angle equilibrium deviation q95/q99/max: "
            f"{quantile(stage03b_angle_eq_deviation, 0.95):.4f}/"
            f"{quantile(stage03b_angle_eq_deviation, 0.99):.4f}/"
            f"{stage03b_angle_eq_deviation.max():.4f} degrees\n\n"
        )

        handle.write(
            "Stage04 remains unauthorized until this "
            "diagnostic is interpreted.\n"
        )

    print(
        "Day021 Stage03b structural-equilibration "
        "diagnostic completed."
    )

    print(
        f"Serious instability signatures: "
        f"{serious}"
    )

    print(
        f"Harmful non-finite values: "
        f"{harmful}"
    )

    print(
        f"Non-finite XVG values: "
        f"{xvg_nonfinite}"
    )

    print(
        "HBN incremental RMS/median/q95/q99/max: "
        f"{stage03b_incremental['rms']:.8f}/"
        f"{stage03b_incremental['median']:.8f}/"
        f"{stage03b_incremental['q95']:.8f}/"
        f"{stage03b_incremental['q99']:.8f}/"
        f"{stage03b_incremental['maximum']:.8f} nm"
    )

    print(
        "HBN aligned RMS/max: "
        f"{stage03b_incremental['aligned_rms']:.8f}/"
        f"{stage03b_incremental['aligned_max']:.8f} nm"
    )

    print(
        "HBN atoms above "
        "0.05/0.08/0.10/0.12/0.15 nm: "
        f"{stage03b_incremental['above_0p05']}/"
        f"{stage03b_incremental['above_0p08']}/"
        f"{stage03b_incremental['above_0p10']}/"
        f"{stage03b_incremental['above_0p12']}/"
        f"{stage03b_incremental['above_0p15']}"
    )

    print(
        f"Stage03/Stage03b top-20 overlap: "
        f"{top20_overlap}/20"
    )

    print(
        "HBN bond length mean/min/max: "
        f"{stage03b_bonds.mean():.8f}/"
        f"{stage03b_bonds.min():.8f}/"
        f"{stage03b_bonds.max():.8f} nm"
    )

    print(
        "HBN bond-equilibrium deviation "
        "q95/q99/max: "
        f"{quantile(stage03b_bond_eq_deviation, 0.95):.8f}/"
        f"{quantile(stage03b_bond_eq_deviation, 0.99):.8f}/"
        f"{stage03b_bond_eq_deviation.max():.8f} nm"
    )

    print(
        "HBN angle-equilibrium deviation "
        "q95/q99/max: "
        f"{quantile(stage03b_angle_eq_deviation, 0.95):.4f}/"
        f"{quantile(stage03b_angle_eq_deviation, 0.99):.4f}/"
        f"{stage03b_angle_eq_deviation.max():.4f} deg"
    )

    print(
        f"Outlier atom names above 0.08 nm: "
        f"{dict(displaced_names)}"
    )

    print(
        f"Structural screen: "
        f"{structural_screen}"
    )

    if review_flags:
        print(
            "Review flags: "
            + " | ".join(
                review_flags
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
        "Stage04 authorized: NO"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )


if __name__ == "__main__":
    main()
